from collections import namedtuple
from typing import TextIO

class REDItor:
    def __init__(self) -> None:
        self.header = [
            "Region", 
            "Position", 
            "Reference", 
            "Strand", 
            "Coverage_q25", 
            "MeanQ", 
            "BaseCount", 
            "AllSubs", 
            "Frequency",
            "gCoverage_q25", 
            "gMeanQ", 
            "gBaseCount", 
            "gAllSubs", 
            "gFrequency"
        ]
        self.field_func = {
            "Region": str, 
            "Position": int, 
            "Reference": str, 
            "Strand": int, 
            "Coverage_q25": int, 
            "MeanQ": float, 
            "BaseCount": lambda x: eval(x), 
            "AllSubs": lambda x: x.split(" "), 
            "Frequency": float,
            "gCoverage_q25": int, 
            "gMeanQ": float, 
            "gBaseCount": lambda x: eval(x), 
            "gAllSubs": lambda x: x.split(" "), 
            "gFrequency": float
        }
        self.REDItuple = namedtuple(
            typename="REDItuple",
            field_names=self.header
        )
        pass
    
    def parse(
        self,
        line: str
        ) -> namedtuple:
        val = []
        for (h, func), raw in zip(
            self.field_func.items(),
            line.strip().split("\t")
        ):
            val.append(func(raw))        
        return self.REDItuple._make(val)
        
    def index(
        self,
        filename: str,
        idx_func: callable = lambda x: x.strip().split("\t")[:2],
        targets: set = None
        ) -> dict:
        with open(filename) as handle:
            idx = {}
            line = handle.readline()
            while line:
                if targets:
                    if tuple(idx_func(line)) in targets:
                        idx[tuple(idx_func(line))] = handle.tell()
                else:
                    idx[tuple(idx_func(line))] = handle.tell()
                line = handle.readline()
        return idx
        
# class REDItable(object):
#     def __init__(
#         self,
#         fname: str,
#         fidx: str = None
#         ) -> None:
#         self.fname = fname
#         self.fidx = fidx if fidx else f"{fname}.idx"
        
#         pass
        
        
if __name__ == "__main__":
    from time import time
    r = REDItor()
    handle = open("pyout/archive-20230104/egr-tsp4-xr.txt")
    line = handle.readline()
    t = r.parser(line)
    time()
    idx = r.index("pyout/archive-20230104/egr-tsp4-xr.txt")
    time()
    print(t)
    print("")

