import pandas as pd
from rnavigate import data, transcriptomics
import requests
from pathlib import Path


def download_eclip_peaks(outpath, assembly="GRCh38"):
    """download eCLIP bed files from ENCODE database

    Args:
        assembly (string, optional): reference genome ("h19" or "GRCh38")
            Defaults to "GRCh38"
        outpath (string, optional): output directory path
            Defaults to "eCLIP_downloads".
    """
    outpath = Path(outpath)
    if not Path.is_dir(outpath):
        Path.mkdir(outpath, parents=True)
        print(f"Created {outpath} folder.")

    print("Searching ENCODE database for eCLIP experiments")
    accession = []
    bed_accession = []

    # Force return from the server in JSON format
    headers = {"accept": "application/json"}
    # GET accession number for experiments
    url_bed = (
        "https://www.encodeproject.org/search/?"
        "type=Experiment&"
        "assay_title=eCLIP&"
        "frame=object&"
        "limit=all&"
        "files.file_type=bed+narrowPeak"
    )
    response = requests.get(url_bed, headers=headers, timeout=21.1)
    # Load json response as a Python dictionary and append accession to list
    experiment_dict = response.json()
    for entry in experiment_dict["@graph"]:
        accession.append(entry["accession"])

    # GET accession number for IDR bed files, store in a dict.
    # Key is experiment acc and values are bed files acc for that experiment.
    for exp in accession:
        url2 = (
            "https://www.encodeproject.org/search/?"
            "type=File&"
            f"dataset=/experiments/{exp}/&"
            "file_format=bed&"
            "format=json&"
            "frame=object&"
            "limit=all"
        )
        response2 = requests.get(url2, headers=headers, timeout=21.1)
        bed_dict = response2.json()
        for i in bed_dict["@graph"]:
            # choose correct assembly and choose replicable peaks
            # if want a specific replicate then change
            # i["biological_replicates"] == "rep_1" or "rep_2"
            if (len(i["biological_replicates"]) == 2) and (i["assembly"] == assembly):
                bed_accession.append(i["accession"])

    # download all IDR bed files in accession
    for acc in bed_accession:
        bedgzfile = requests.get(
            f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.bed.gz",
            timeout=21.1,
        )
        with open(outpath / f"{acc}.bed.gz", "wb") as outfile:
            outfile.write(bedgzfile.content)
    print(f"Download finished. All files are stored in {outpath} dir.")


def create_eclip_table(inpath, outpath):
    # create a table to look up eCLIP files using target and cell type.
    codes = {"accession": [], "target": [], "cell_line": []}
    for file in Path(inpath).iterdir():
        if not file.match("*.bed.gz"):
            continue
        df = pd.read_table(file, compression="gzip", header=None)
        label = df.iloc[0, 3]
        label = label.split("_")
        if len(label) == 1:
            label = ["", "", ""]
        elif len(label) == 4:
            label = [label[0], f"{label[1]}_{label[2]}", label[3]]
        codes["accession"].append(file.name.rstrip(".bed.gz"))
        codes["target"].append(label[0])
        codes["cell_line"].append(label[1])
    codes = pd.DataFrame(codes)
    k562 = codes.loc[codes["cell_line"] == "K562", ["target", "accession"]]
    k562.columns = ["target", "K562"]
    hepg2 = codes.loc[codes["cell_line"] == "HepG2", ["target", "accession"]]
    hepg2.columns = ["target", "HepG2"]
    df = pd.merge(k562, hepg2, how="outer", on="target")
    df[["target", "K562", "HepG2"]].to_csv(
        Path(outpath) / "eclip_codes.txt", index=False, sep="\t"
    )


class eCLIPDatabase:  # pylint disable=invalid-name
    def __init__(self, inpath):
        self.path = Path(inpath)
        self.eclip_codes = pd.read_table(self.path / "eclip_codes.txt")
        self.eclip_data = self.get_eclip_data()

    def get_cell_target_data(self, cell_line, target):
        try:
            cell_target = self.eclip_data[cell_line]
        except KeyError:
            raise ValueError(f"Cell line {cell_line} not found.")
        try:
            cell_target = cell_target[target]
        except KeyError:
            raise ValueError(f"Target {target} not found for cell line {cell_line}.")
        return cell_target

    def get_eclip_data(self):
        eclip_data = {"HepG2": {}, "K562": {}}
        hepg2_rows = ~self.eclip_codes["HepG2"].isnull()
        k562_rows = ~self.eclip_codes["K562"].isnull()
        for _, row in self.eclip_codes[hepg2_rows].iterrows():
            eclip_data["HepG2"][row["target"]] = transcriptomics.NarrowPeak(
                str(self.path / f"{row['HepG2']}.bed.gz")
            )
        for _, row in self.eclip_codes[k562_rows].iterrows():
            eclip_data["K562"][row["target"]] = transcriptomics.NarrowPeak(
                str(self.path / f"{row['K562']}.bed.gz")
            )
        return eclip_data

    def get_eclip_density(self, transcript, cell_line, targets=None):
        eclip_data = self.eclip_data[cell_line]
        if targets is not None:
            eclip_data = {target: eclip_data[target] for target in targets}
        eclip_density = data.Profile(
            input_data=transcript.coordinate_df.eval("eCLIP_density = 0"),
            metric="eCLIP_density",
        )
        for value in eclip_data.values():
            values = value.get_profile(transcript).data.eval("value != 0.0")
            eclip_density.data["eCLIP_density"] += values
        return eclip_density

    def get_annotation(self, transcript, cell_line, target, **kwargs):
        cell_target = self.get_cell_target_data(cell_line, target)
        return cell_target.get_annotation(transcript, **kwargs)

    def get_profile(self, transcript, cell_line, target):
        cell_target = self.get_cell_target_data(cell_line, target)
        return cell_target.get_profile(transcript)

    def print_all_peaks(self, transcript):
        for target in self.eclip_codes["target"]:
            for cell_line in ["K562", "HepG2"]:
                self.print_peaks(transcript, cell_line, target)

    def print_peaks(self, transcript, cell_line, target):
        try:
            annotation = self.get_annotation(
                transcript=transcript, cell_line=cell_line, target=target
            )
            if len(annotation.data) != 0:
                df = annotation.data.sort_values("start")
                spans = [f"{mn}-{mx}" for _, (mn, mx) in df.iterrows()]
                print(f"{cell_line:<8} {target:<8} {' '.join(spans)}")
        except ValueError as e:
            pass
