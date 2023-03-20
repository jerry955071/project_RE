# IMPORT
# ==============================================================================
from pathlib import Path
from collections import deque
import sys
import os


# VARS
# ==============================================================================
FIN = sys.argv[1]
PATH_OUT = Path(sys.argv[2])
WORKDIR = Path("/home/b05b01002/HDD3/project_RE")

# CHDIR
# ==============================================================================
os.chdir(WORKDIR)


# 1. REDItools --(clean)--> Archive
# ==============================================================================
# get files paths
fin = Path(FIN)
fout = PATH_OUT / fin.name

# create output directory
os.makedirs(fout.parent, exist_ok=True) 

# open files
handle_in = open(fin, "r")
handle_out = open(fout, "w")

# print start
print(f"Processing file: {str(fin)}", file=sys.stderr)

# clean lines
d = deque(["", ""], maxlen=2)
n = 0
for line in handle_in:
    n += 1
    if n % 10000000 == 0:
        print(f"Read {n} lines from '{str(fin)}'", file=sys.stdout)
    d.append(line)
    if d[0] == "":
        continue
    else:
        chr0, pos0 = d[0].split("\t")[:2]
        chr1, pos1 = d[1].split("\t")[:2]
        if (chr0 != chr1) or (pos0 != pos1):
            if d[0].split("\t")[9] != "-":
                foo = handle_out.write(d[0])
        if chr0 != chr1:
            print(f"Processing region: {chr1}", file=sys.stdout)

# write the last line in deque
foo = handle_out.write(d[1])
    
# close file handles
handle_in.close()
handle_out.close()

# print finished
print(f"Finished processing file: {str(fin)}", file=sys.stderr)
