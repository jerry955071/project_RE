from mods.Location import Location
from mods.REDItor import REDItor
from typing import TextIO

def skimming(
        handle: TextIO,
        target: Location,
        ref_order: list
    ):
    while True:
        # try read line
        try:
            streampos = handle.tell()
            newline = handle.readline()
        except StopIteration:
            return None
        
        # check position
        t = REDItor().parser(newline)
        new_loc = Location(t.Region, t.Position, name="<name>", order=ref_order)
        if new_loc == target:
            return t
        if new_loc > target:
            handle.seek(streampos)
            return None
        
        
if __name__ == "__main__":
    ref_order = [
        line.split("\t")[0] 
        for line in open("genomic_data/Egr/Egr.fa.fai", "r")
    ]
    Location.set_ref_order = ref_order
    r = REDItor()
    handle = open("pyout/archive-20230104/egr-tsp4-xr.txt")
    skimming(handle, Location("Chr01", 20000, ref_order), ref_order)
    
    print("")
