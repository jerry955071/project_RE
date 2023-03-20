# conda: env/re.yaml
# %%
from multiprocessing import Pool
from pathlib import Path
import pandas as pd
import json
import sys
import os

# %%
# get run paramters
# params = json.load(open("configs/differential_editing/Ptr.json"))
params = json.load(open(sys.argv[1]))
os.makedirs(params["tmpdir1"], exist_ok=True)

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


def main(params, reditables, target_positions):
    for idx, (reg, pos) in target_positions.iterrows():
        allsubs = set()
        basecounts = []
        rname = []
        for sample, df in reditables.items():
            # try get row
            try:
                # record observed substitution types
                if pd.isna(df.loc[(reg, pos), "AllSubs"])[0]:
                    pass
                else:
                    for i in df.loc[(reg, pos), "AllSubs"][0].split(" "):
                        allsubs.add(i)
                
                # record basecounts
                rname.append(sample)
                basecounts.append(
                    df.loc[(reg, pos), "BaseCount"][0] + \
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
            fout = f"{params['tmpdir1']}/{reg}_{pos}_{subtype}.csv"
            pd.DataFrame(
                basecounts,
                index=rname,
                columns=["A","C","G","T", "Tissue"]
            ).to_csv(fout)

# %%
# partitioning target_positions
nrow = target_positions.shape[0]
intervals = [i for i in range(0, nrow, nrow // params["threads"])]
intervals.append(nrow)
partitions = []
for i in range(params["threads"]):
    partitions.append(target_positions.iloc[intervals[i]:intervals[i+1],:])


# multiprocessing on each partition of the target positions
# %%
if __name__ == "__main__":
    with Pool(processes=params["threads"]) as p:
        p.starmap(
            main, 
            [(params, reditables, nth_part) for nth_part in partitions]
        )
    

# Timing:
# Running on 345950 target positions with 57 samples using 1 thread 
# 56:27.92 total

# Running on 345950 target positions with 57 samples using 20 processes
# 5710.05s user 333.63s system 731% cpu 13:46.00 total
# %%
