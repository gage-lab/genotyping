import requests, json, time, secrets, re, logging, sys
from pathlib import Path
from datetime import datetime
from urllib.request import urlretrieve
from urllib.parse import urlparse
from enum import IntEnum
from zipfile import ZipFile


# Job state constants matching server API values
class JobState(IntEnum):
    QUEUED = 1
    RUNNING = 2
    EXPORTING = 3
    SUCCESS = 4
    FAILED = 5
    CANCELLED = 6
    RETIRED = 7
    RUNNING_PHASE2 = 8
    FAILED_PHASE2 = 9
    DELETED = 10

def setup_logger(name="imputation", log_file=None, level="INFO"):
    """Setup logger with both console and optional file output."""
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper()))
    logger.handlers.clear()

    # Console handler with timestamp formatting
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(
        logging.Formatter(
            "[%(asctime)s] [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S"
        )
    )
    logger.addHandler(console)

    # Optional file handler (plain text, no colors)
    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(
            logging.Formatter(
                "[%(asctime)s] [%(levelname)-8s] %(message)s", "%Y-%m-%d %H:%M:%S"
            )
        )
        logger.addHandler(file_handler)
        logger.info(f"Logging to file: {log_file}")

    return logger


class ImputationClient:
    """Client for interacting with Michigan/NIH imputation servers."""

    # Server configuration: (base_url, submit_endpoint)
    SERVERS = {
        "michigan": (
            "https://imputationserver.sph.umich.edu/api/v2",
            "/jobs/submit/imputationserver2",
        ),
        # TODO: test NIH server
        # "nih": (
        #     "https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2",
        #     "/jobs/submit/minimac4",
        # ),
    }

    def __init__(self, server: str, token: str, logger=None):
        """Initialize client with server config and authentication."""
        self.logger = logger or logging.getLogger(__name__)
        self.url, self.submit_endpoint = self.SERVERS[server.lower()]

        # Use session for connection reuse and persistent headers
        self.session = requests.Session()
        self.session.headers.update({"X-Auth-Token": token})
        self.logger.info(f"Initialized client for {server} server: {self.url}")

    def _request(self, method: str, endpoint: str, **kwargs):
        """Make HTTP request with consistent error handling and logging."""
        self.logger.debug(f"Making {method} request to {endpoint}")
        r = self.session.request(method, f"{self.url}{endpoint}", **kwargs)

        # Handle common API errors
        if r.status_code == 401:
            self.logger.error("Authentication failed - bad or expired API token")
            raise ValueError("Bad API token")
        if r.status_code == 404:
            self.logger.error(f"API endpoint not found: {endpoint}")
            raise ValueError("Invalid API URL")

        r.raise_for_status()  # Raise for other HTTP errors
        self.logger.debug(f"API response: {r.status_code}")
        return r

    def get_jobs(self):
        """Fetch all jobs from the imputation server."""
        self.logger.debug("Fetching all jobs from server")
        jobs = self._request("GET", "/jobs").json()["data"]
        self.logger.debug(f"Retrieved {len(jobs)} jobs from server")
        return jobs

    def get_job_info(self, job_id: str):
        """Get detailed information for a specific job."""
        self.logger.debug(f"Fetching info for job {job_id}")
        return self._request("GET", f"/jobs/{job_id}").json()

    def wait_for_queue_slot(self, max_jobs=2):
        """Wait until there are available slots to submit new jobs."""
        self.logger.info(
            f"Checking for available queue slots (max {max_jobs} concurrent jobs)"
        )

        while True:
            # Get jobs that haven't completed yet
            incomplete = [j for j in self.get_jobs() if j["state"] < JobState.SUCCESS]

            # Check if we can submit (within concurrent limit)
            if len(incomplete) <= max_jobs:
                msg = (
                    f"✓ Ready to submit. Currently {len(incomplete)} job(s) in progress"
                    if incomplete
                    else "✓ No jobs currently running. Ready to submit"
                )
                self.logger.info(msg)
                break

            self.logger.warning(
                f"⏳ {len(incomplete)} jobs are queued or running. Waiting..."
            )

            # Show detailed queue status
            running = [j for j in incomplete if j["state"] > JobState.QUEUED]
            if not running:  # Only queued jobs
                min_pos = min(j["positionInQueue"] for j in incomplete)
                self.logger.info(f"📊 Lowest queue position: {min_pos}")
            else:  # Some jobs running
                self.logger.info(f"🏃 {len(running)} jobs currently running")

            self.logger.info("Waiting 10 minutes before checking again...")
            time.sleep(600)  # Wait 10 minutes between checks

    def submit_job(self, vcf_files: list, build: str, refpanel: str, population: str) -> dict:
        """Submit a new imputation job with the provided VCF files."""
        message = f"Preparing to submit job with {len(vcf_files)} VCF file(s):"
        for vcf in vcf_files:
            message += f"\n  - {vcf}"
        message += "\n\nJob settings:"
        message += f"\n  - Build: {build}"
        message += f"\n  - Reference panel: {refpanel}"
        message += f"\n  - Population: {population}"
        self.logger.info(message)

        # Prepare job configuration with secure password and timestamp
        data = {
            "password": secrets.token_urlsafe(36),  # Generate secure password
            "job-name": f"job_{datetime.now():%Y-%m-%d.%H%M}",
            "build": build,
            "refpanel": refpanel,
            "population": population,
        }
        self.logger.info(f"Job name: {data['job-name']}")

        # Submit with proper file handling
        with open_files(vcf_files) as files:
            self.logger.info("🚀 Submitting job to server...")
            r = self._request("POST", self.submit_endpoint, files=files, data=data)

        # Enhance response with our local data
        result = r.json()
        result.update(
            {
                "password": data["password"],
                "job-name": data["job-name"],
                "settings": {
                    k: v for k, v in data.items() if k not in ["password", "job-name"]
                },
            }
        )
        self.logger.info(f"✅ Job submitted successfully! Job ID: {result['id']}")
        return result

    def monitor_job(self, job_id: str) -> dict:
        """Monitor job progress until completion."""
        self.logger.info(f"👀 Monitoring job {job_id}...")
        info = self.get_job_info(job_id)
        submission_time = datetime.fromtimestamp(info["submittedOn"] / 1000)
        self.logger.info(
            f"📍 Initial queue position: {info.get('positionInQueue', 'N/A')}"
        )

        _fmt_delta = lambda x: str(datetime.now() - x).split(".")[0]

        # Monitor until job completes
        while info["state"] < JobState.SUCCESS:
            # Format status messages based on current state
            status_msgs = {
                JobState.QUEUED: f"⏳ {job_id} waiting for {_fmt_delta(submission_time)} (Queue position: {info.get('positionInQueue', 'Unknown')})",
                JobState.RUNNING: f"🏃 {job_id} running for {_fmt_delta(datetime.fromtimestamp(info['startTime'] / 1000))}",
                JobState.EXPORTING: f"📦 {job_id} exporting results",
            }

            if info["state"] in status_msgs:
                self.logger.info(status_msgs[info["state"]])
            else:
                self.logger.error(f"💀 {job_id} died with state {info['state']}")
                raise Exception(f"{job_id} died with state {info['state']}")

            # Different polling intervals based on job state
            time.sleep(60 if info["state"] == JobState.QUEUED else 20)
            info = self.get_job_info(job_id)

        self.logger.info(f"🎉 Job {job_id} completed with state {info['state']}")
        return info

    def download_results(self, job_info: dict, output_path: str, password: str):
        """Download all result files from completed job."""
        # Check if results are available
        if job_info["state"] in [JobState.RETIRED, JobState.DELETED]:
            self.logger.warning("⚠️  Job results are not available for download")
            return

        # Prepare download URLs and file list
        friendly_name = re.sub(
            r"_submitted20\d\d-\d\d-\d\d\.\d+$", "", job_info["name"]
        )
        base_url = (
            f"https://{urlparse(self.url).netloc}/share/results/{{hash}}/{{name}}"
        )

        # Calculate total files for logging
        total_files = sum(len(param["files"]) for param in job_info["outputParams"] if "files" in param)
        self.logger.info(
            f"📁 Prepared {total_files} files for download across {len([p for p in job_info['outputParams'] if 'files' in p])} categories"
        )

        # Download each category of files directly
        for param in job_info["outputParams"]:
            if "files" not in param:
                continue
            
            desc = param["description"]
            self.logger.info(f"📥 Downloading {desc} for {friendly_name}")
            
            for file_data in param["files"]:
                dest = Path(output_path) / file_data["name"]
                # Ensure output directory exists
                dest.parent.mkdir(parents=True, exist_ok=True)

                self.logger.info(
                    f"⬇️  Starting download: {file_data['name']} ({file_data['size']})"
                )
                
                # Download file
                download_url = base_url.format(**file_data)
                urlretrieve(download_url, dest)

                if file_data["name"].endswith(".zip"):
                    self.logger.info(f"🔓 Unzipping {file_data['name']}")
                    with ZipFile(dest) as zip_ref:
                        zip_ref.extractall(dest.parent, pwd=password.encode("utf-8"))
                    dest.unlink()
                
                self.logger.info(f"✅ Finished downloading {file_data['name']}")

def open_files(vcf_files):
    """Context manager for safe file handling - ensures files are always closed."""
    from contextlib import contextmanager

    @contextmanager
    def _open():
        files = [open(vcf, "rb") for vcf in vcf_files]
        try:
            yield [("files", f) for f in files]  # Format for requests
        finally:
            [f.close() for f in files]  # Always close files

    return _open()


def handle_completion(job_info: dict, job_id: str, logger):
    """Handle job completion status and exit appropriately."""
    state = job_info["state"]

    # Success states
    if state in [JobState.SUCCESS, JobState.RUNNING_PHASE2]:
        logger.info(f"🎉 {job_id} completed successfully!")
        return

    # Failure states with appropriate messages
    failures = {
        JobState.FAILED: "failed. Check logs",
        JobState.FAILED_PHASE2: "failed. Check logs",
        JobState.CANCELLED: "cancelled. Please rerun",
        JobState.RETIRED: "retired. Please rerun",
        JobState.DELETED: "deleted. Please rerun",
    }

    if state in failures:
        logger.error(f"❌ {job_id} was {failures[state]}")
        sys.exit(1)

    # Unknown state
    logger.error(f"❓ {job_id} status not recognized. Status = {state}")
    sys.exit(1)


if __name__ == "__main__":
    # Setup logging with file output from Snakemake
    logger = setup_logger("imputation", log_file=str(snakemake.log), level="INFO")

    # Log workflow startup information
    logger.info("🚀 Starting imputation workflow")
    output_dir = Path(snakemake.output.qc_report).parent
    logger.info(f"Output directory: {output_dir}")

    try:
        # Initialize client and wait for available queue slot
        client = ImputationClient(
            snakemake.params["server"], snakemake.params["token"], logger
        )
        client.wait_for_queue_slot()

        # Submit job and log submission details
        submission = client.submit_job(
            vcf_files=snakemake.input.vcf,
            build=snakemake.params["build"],
            refpanel=snakemake.params["refpanel"],
            population=snakemake.params["population"],
        )
        logger.info(f"📋 {submission['message']}")
        logger.info(f"🆔 Job ID: {submission['id']}")

        # Monitor job until completion
        final_info = client.monitor_job(submission["id"])

        # Download results to output directory
        client.download_results(final_info, str(output_dir), password=submission["password"])

        # Save complete job information for reference
        job_info_path = Path(snakemake.output.job_json)
        job_info_path.write_text(json.dumps(final_info, indent=2))
        logger.info(f"💾 Saved job information to {job_info_path}")
    
        # Handle final job status and exit appropriately
        handle_completion(final_info, submission["id"], logger)
        logger.info("🏁 Workflow completed successfully!")

    except Exception as e:
        # Log full error details and exit with error code
        logger.error(f"💥 Workflow failed: {str(e)}", exc_info=True)
        sys.exit(1)
