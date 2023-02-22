# conda: env/re.yaml
# %%
from pathlib import Path
import pandas as pd
import json
import sys
import os

# %%
# get run paramters
params = json.load(open("configs/differential_editing/Ptr.json"))
os.makedirs(params["tmpdir"])
# params = json.load(open(sys.argv[1]))

# %%
# import dataframes
samples = {}
reditables = {}
for tissue, sample_paths in params["samples"].items():
    samples[tissue] = []
    for sample_path in sample_paths:
        sample_path = Path(sample_path)
        reditables[sample_path.name] = pd.read_csv(
            sample_path,
            sep="\t",
            header=None,
            names=list(params["redithead"].keys()),
            index_col=["Region", "Position"],
            usecols=[
                "Region",
                "Position",
                "Reference",
                "Strand",
                "BaseCount",
                "AllSubs"
            ],
            converters={"BaseCount":eval},
            dtype=params["redithead"],
            na_values="-"
        )
        samples[tissue].append(str(sample_path.name))

# %%
target_positions = pd.read_csv(
    params["target_positions"],
    sep="\t",
    header=None,
    names=("Region", "Position")
)

# %%
N = 2000
M = 0
for idx, (reg, pos) in target_positions.iterrows():
    if M > N:
        break
    else:
        M += 1
    allsubs = set()
    basecounts = []
    rname = []
    for sample, df in reditables.items():
        # try get row
        try:
            # record observed substitution types
            if pd.isna(df.loc[(reg, pos), "AllSubs"]):
                pass
            else:
                for i in df.loc[(reg, pos), "AllSubs"].split(" "):
                    allsubs.add(i)
            
            # record basecounts
            rname.append(sample)
            basecounts.append(
                df.loc[(reg, pos), "BaseCount"] + \
                [
                    {
                        "x": "xylem",
                        "p": "phloem",
                        "l": "leaf",
                        "s": "shoot"
                    }[sample[-6]]
                ]
            )
            
        except KeyError:
            pass
    
    for subtype in allsubs:
        fout = f"{params['tmpdir']}/{reg}_{pos}_{subtype}.csv"
        pd.DataFrame(
            basecounts,
            index=rname,
            columns=["A","C","G","T", "Tissue"]
        ).to_csv(fout)

# %%
