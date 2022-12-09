# IMPORT
# ==============================================================================
from pathlib import Path
from collections import deque, defaultdict
from typing import TextIO
from mods.Location import Location
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
        target: Location,
        order: list
    ):
    # read one line
    newline = next(handle)      
    while True:
        # if EOF
        if (newline == None) or (newline == ""):
            return None, None
        
        # check position
        regq, posq = newline.split("\t")[:2]
        posq = int(posq)
        new_loc = Location(regq, posq, order)
        if new_loc >= target:
            return newline, new_loc
        
        # else repeat
        newline = next(handle)


# ==============================================================================
# While
ref_order = [line.split("\t")[0] for line in open(FAI[_SPECIES], "r")]
fin_lst = os.listdir(ARCHIVE_PATH[_SPECIES])
fin_handles = {f: open(ARCHIVE_PATH[_SPECIES] / f) for f in fin_lst}
fout_handles = {f: open(HEXPR_PATH[_SPECIES] / f, "w") for f in fin_lst}
fout_handles = {f: open(HEXPR_PATH[_SPECIES] / f, "a") for f in fin_lst}
line_buffer = {f: tuple([None, Location("Chr01", -1, ref_order)]) for f in fin_lst}

_EOF = False
target_loc = Location("Chr01", 1, ref_order)
while not _EOF:
    # scan for target location reg:pos
    print(f"[Region] {target_loc}", file=sys.stderr, end="")    
    
    # for each file in fin_lst
    _TARGET_UNCHANGED = True
    for file in fin_lst:
        # update line if line.loc < target.loc
        if line_buffer[file][1] < target_loc:
            newline, new_loc = skimming(fin_handles[file], target_loc, ref_order)
            line_buffer[file] = (newline, new_loc)

        # if any file ends
        if newline == None:
            _EOF = True
            break

        # update target location if newline passes target location
        if new_loc > target_loc:
            print("\tSkipped", file=sys.stderr)
            _TARGET_UNCHANGED = False
            target_loc = new_loc
            break
    
    # if target unchanged
    if _TARGET_UNCHANGED:
        print("\tMatched", file=sys.stderr)
        for k, v in line_buffer.items():
            fout_handles[k].write(v[0])
        target_loc += 1

for k, v in fin_handles.items():
    v.close() 
for k, v in fout_handles.items():
    v.close() 