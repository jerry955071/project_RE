# Created: 2022-12-07  
# ==============================================================================
# 1. REDItools output --clean--> ptr-edsites-20221207 
# 2. Identify tissue-specific RNA editing events  
# ==============================================================================


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


# # GET FILES
# # ==============================================================================
# Q = Query()
# mydb = TinyDB(MYDB_PATH / "db.json")
# for spe in _SPECIES:
#     records = mydb.search((Q.species == spe) & (Q.type == "RNA"))
#     _REDIOUT += [
#             "processed/reditools2/%s/Final/%s/%s.txt" % (
#             r["species"],
#             r["tissue"],
#             r["id"]
#             ) for r in records
#         ]


# # NOTE: DONE
# # 1. REDItools --(clean)--> Archive
# # ==============================================================================
# for f in _REDIOUT:
#     # open files
#     fin = Path(f)
#     fout = ARCHIVE_PATH / fin.name
#     os.makedirs(fout.parent, exist_ok=True) 
#     handle_in = open(f, "r")
#     handle_out = open(fout, "w")
#     print(f"Processing file: {str(fin)}", file=sys.stderr)

#     # clean lines
#     d = deque(["", ""], maxlen=2)
#     n = 0
#     for line in handle_in:
#         n += 1
#         if n % 10000 == 0:
#             print(f"Read {n} lines from '{str(fin)}'", file=sys.stdout)
#         d.append(line)
#         if d[0] == "":
#             continue
#         else:
#             chr0, pos0 = d[0].split("\t")[:2]
#             chr1, pos1 = d[1].split("\t")[:2]
#             if (chr0 != chr1) or (pos0 != pos1):
#                 foo = handle_out.write(d[0])
#             if chr0 != chr1:
#                 print(f"Processing region: {chr1}", file=sys.stdout)
#     foo = handle_out.write(d[1])
    
#     # end of loop
#     handle_in.close()
#     handle_out.close()
    

# 2. Identify tissue-specific RNA editing events  
# ==============================================================================
# List files
os.listdir(ARCHIVE_PATH)

# `skimming`
def skimming(
        handle: TextIO,
        region: str,
        pos: int
    ):
    newline = next(handle)      
    while True:
        if not newline:
            return None, None, None
        
        regq, posq = newline.split("\t")[:2]
        posq = int(posq)
        
        if (regq == region) and (posq >= pos):
            return newline, regq, posq
        
        if regq != region:
            return newline, regq, posq
        
        newline = next(handle)

# base stage
PATH_NEW = Path("hexpr/") / str(round(time.time()))
os.makedirs(PATH_NEW, exist_ok=True)
fin_handles = {f: open(ARCHIVE_PATH / f) for f in os.listdir(ARCHIVE_PATH)}
fout_handles = {f: open(PATH_NEW / f, "w") for f in os.listdir(ARCHIVE_PATH)}
fout_handles = {f: open(PATH_NEW / f, "a") for f in os.listdir(ARCHIVE_PATH)}
line_buffer = {f: tuple([None, None, None]) for f in os.listdir(ARCHIVE_PATH)}

_EOF = False
reg = "Chr01"
pos = 1
while not _EOF:
    print(f"[Region] {reg}:{pos}", file=sys.stderr, end="")
    # if any file ends
    if _EOF:
        break
    
    _TARGET_UNCHANGED = True
    for file in fin_handles.keys():
        if not ((line_buffer[file][1] == reg) and (line_buffer[file][2] == pos)):
            # skimming
            newline, new_reg, new_pos = skimming(fin_handles[file], reg, pos)
            line_buffer[file] = (newline, new_reg, new_pos)

        # if any file ends
        if not newline:
            _EOF = True
            break

        # update target location if newline passes target location
        if not ((new_reg == reg) and (new_pos == pos)):
            print("\tSkipped", file=sys.stderr)
            _TARGET_UNCHANGED = False
            reg = new_reg
            pos = new_pos
            break
    
    # if target unchanged
    if _TARGET_UNCHANGED:
        print("\tMatched", file=sys.stderr)
        for k, v in line_buffer.items():
            fout_handles[k].write(v[0])
        pos += 1

for k, v in fin_handles.items():
    v.close() 
for k, v in fout_handles.items():
    v.close() 