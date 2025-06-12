from pathlib import Path
import pandas as pd

# corresponding supporting files for each array type
ARRAY_CONFIG = {
    "GSA-24v3-0_A2": {
        "bpm_manifest_file": "resources/illumina_array_support/GSA-24v3-0_A2.bpm",
        "egt_cluster_file": "resources/illumina_array_support/GSA-24v3-0_A1_ClusterFile.egt",
        "csv_manifest_file": "resources/illumina_array_support/GSA-24v3-0_A2.csv",
        "ref": "resources/Homo_sapiens_assembly38.fasta",
    },
    "InfiniumCoreExome-24v1-4_A1": {
        "bpm_manifest_file": "resources/illumina_array_support/InfiniumCoreExome-24v1-4_A1.bpm",
        "egt_cluster_file": "resources/illumina_array_support/InfiniumCoreExome-24v1-4_A1_ClusterFile.egt",
        "csv_manifest_file": "resources/illumina_array_support/InfiniumCoreExome-24v1-4_A1.csv",
        "ref": "resources/Homo_sapiens_assembly19.fasta",
    },
}


def read_samplesheet(igm_dir: str) -> tuple[pd.DataFrame, dict]:
    """
    Read the samplesheet from the igm_dir.
    Returns a tuple of a pandas DataFrame and a dictionary of array configuration.
    """
    
    search = [p for p in Path(igm_dir).rglob("*.csv") if not p.name.startswith(".")]
    if len(search) == 0:
        raise FileNotFoundError(f"No samplesheet found in {igm_dir}")
    samplesheet = search[0]

    print(f"Reading samplesheet from {samplesheet}")
    with open(samplesheet, "r") as f:
        for i, line in enumerate(f):
            if i == 17:
                manifest = line.split(",")[1]
                break

    if manifest not in ARRAY_CONFIG:
        raise ValueError(f"Manifest {manifest} not found in ARRAY_CONFIG")

    array_config = ARRAY_CONFIG[manifest]
    samples = pd.read_csv(
        samplesheet,
        dtype={"SentrixBarcode_A": str, "SentrixPosition_A": str},
        skiprows=22,
    )
    samples["sentrix_id"] = (
        samples["SentrixBarcode_A"] + "_" + samples["SentrixPosition_A"]
    )

    # check that manifest matches bpm_manifest_file and csv_manifest_file
    if manifest != Path(array_config["bpm_manifest_file"]).stem:
        raise ValueError(
            f"Manifest {manifest} does not match bpm_manifest_file {array_config['bpm_manifest_file']}"
        )
    if manifest != Path(array_config["csv_manifest_file"]).stem:
        raise ValueError(
            f"Manifest {manifest} does not match csv_manifest_file {array_config['csv_manifest_file']}"
        )

    # find idats
    array_config["Grn.idat"] = []
    array_config["Red.idat"] = []
    for _, row in samples.iterrows():
        prefix = Path(igm_dir) / row["SentrixBarcode_A"] / row["sentrix_id"]
        for suffix in ["Grn.idat", "Red.idat"]:
            full_name = str(prefix) + "_" + suffix
            assert Path(full_name).exists(), f"File {full_name} does not exist"
            array_config[suffix].append(full_name)

    # check that idats are found
    if len(array_config["Grn.idat"]) == 0 or len(array_config["Red.idat"]) == 0:
        raise ValueError(f"No idats found for {manifest}")
    elif len(array_config["Grn.idat"]) != len(array_config["Red.idat"]):
        raise ValueError(
            f"Number of green idats ({len(array_config['Grn.idat'])}) and red idats ({len(array_config['Red.idat'])}) do not match"
        )
    elif len(array_config["Grn.idat"]) != len(samples):
        raise ValueError(
            f"Number of idats ({len(array_config['Grn.idat'])}) does not match number of samples ({len(samples)})"
        )

    print(
        f"Found {len(array_config['Grn.idat'])} green idats and {len(array_config['Red.idat'])} red idats"
    )

    # return samples, array_config
    return samples, array_config
