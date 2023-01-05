from pathlib import Path
import sys

#
FIN = Path(sys.argv[1])
PLOC_UNION =  Path(sys.argv[2])
POUT = Path(sys.argv[3])
FOUT = POUT / FIN.name

#
LOC_UNION = set()
for line in open(PLOC_UNION):
    LOC_UNION.add(tuple(line.strip().split("\t")))
    
#
open(FOUT, "w").close()
with (
  open(FIN, "r") as handle_in,
  open(FOUT, "a") as handle_out  
):
    for line in handle_in:
        rname, pos = line.strip().split("\t")[:2]
        if (rname, pos) in LOC_UNION:
            handle_out.write(line)