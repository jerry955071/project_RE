# IMPORT
# ==============================================================================
from tinydb import TinyDB, Query
from pathlib import Path
from collections import deque, defaultdict
from typing import TextIO
import time
import sys
import os



# VARS
# ==============================================================================
WORKDIR = Path("/home/b05b01002/HDD3/project_RE")
ARCHIVE_PATH = Path("pyout/archive-20221207")
MYDB_PATH = Path("mydb")
FAI = {
    "Ptr": Path("genomic_data/Ptr/Ptr.fa.fai"),
    "Egr": Path("genomic_data/Egr/Egr.fa.fai")
}
_SPECIES = ["Egr"]
_REDIOUT = []


# CHDIR
# ==============================================================================
os.chdir(WORKDIR)


# GET FILES
# ==============================================================================
Q = Query()
mydb = TinyDB(MYDB_PATH / "db.json")
for spe in _SPECIES:
    records = mydb.search((Q.species == spe) & (Q.type == "RNA"))
    _REDIOUT += [
            "processed/reditools2/%s/Final/%s/%s.txt" % (
            r["species"],
            r["tissue"],
            r["id"]
            ) for r in records
        ]


# NOTE: DONE
# 1. REDItools --(clean)--> Archive
# ==============================================================================
for f in _REDIOUT:
    # open files
    fin = Path(f)
    fout = ARCHIVE_PATH / fin.name
    os.makedirs(fout.parent, exist_ok=True) 
    handle_in = open(f, "r")
    handle_out = open(fout, "w")
    print(f"Processing file: {str(fin)}", file=sys.stderr)

    # clean lines
    d = deque(["", ""], maxlen=2)
    n = 0
    for line in handle_in:
        n += 1
        if n % 10000 == 0:
            print(f"Read {n} lines from '{str(fin)}'", file=sys.stdout)
        d.append(line)
        if d[0] == "":
            continue
        else:
            chr0, pos0 = d[0].split("\t")[:2]
            chr1, pos1 = d[1].split("\t")[:2]
            if (chr0 != chr1) or (pos0 != pos1):
                foo = handle_out.write(d[0])
            if chr0 != chr1:
                print(f"Processing region: {chr1}", file=sys.stdout)
    foo = handle_out.write(d[1])
    
    # end of loop
    handle_in.close()
    handle_out.close()