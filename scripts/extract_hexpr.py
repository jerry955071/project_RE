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
_SPECIES = "Ptr"
WORKDIR = Path("/home/b05b01002/HDD3/project_RE")
MYDB_PATH = Path("mydb")
FAI = {
    "Ptr": Path("genomic_data/Ptr/Ptr.fa.fai"),
    "Egr": Path("genomic_data/Egr/Egr.fa.fai")
}
HEXPR_PATH = {
    "Ptr": Path("pyout/hexpr/Ptr") / str(round(time.time())),
    "Egr": Path("pyout/hexpr/Egr") / str(round(time.time()))
}
ARCHIVE_PATH = {
    "Ptr": Path("pyout/archive-20221207/Ptr"),
    "Egr": Path("pyout/archive-20221207/Egr")

}


# CHDIR
# ==============================================================================
os.chdir(WORKDIR)
os.makedirs(HEXPR_PATH[_SPECIES], exist_ok=True)

# Extract highly expressed regions  
# ==============================================================================
# function `skimming`
def skimming(
        handle: TextIO,
        region: str,
        pos: int
    ):
    # read one line
    newline = next(handle)      
    while True:
        # if EOF
        if not newline:
            return None, None, None
        
        # check position
        regq, posq = newline.split("\t")[:2]
        posq = int(posq)
        
        # same chromosome, position >= targeted position
        if (regq == region) and (posq >= pos):
            return newline, regq, posq
        
        # next chromosome
        if regq != region:
            return newline, regq, posq
        
        newline = next(handle)

# ==============================================================================
# While
fin_lst = os.listdir(ARCHIVE_PATH[_SPECIES])
fin_handles = {f: open(ARCHIVE_PATH[_SPECIES] / f) for f in fin_lst}
fout_handles = {f: open(HEXPR_PATH[_SPECIES] / f, "w") for f in fin_lst}
fout_handles = {f: open(HEXPR_PATH[_SPECIES] / f, "a") for f in fin_lst}
line_buffer = {f: tuple([None, None, None]) for f in fin_lst}

_EOF = False
reg = "Chr01"
pos = 1

while not _EOF:
    print(f"[Region] {reg}:{pos}", file=sys.stderr, end="")    
    
    _TARGET_UNCHANGED = True
    for file in fin_lst:
        if not ((line_buffer[file][1] == reg) and (line_buffer[file][2] >= pos)):
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