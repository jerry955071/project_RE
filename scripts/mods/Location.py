from typing import List

class Location(object):
    def __init__(
            self,
            reg: str,
            pos: int,
            order: List,
            name: str = "<name>"
        ) -> None:
        self.reg = reg
        self.pos = pos
        self.order = order
        self.name = name
        self.val = self.order.index(self.reg) + self.pos * (10 ** -10)
        pass
    
    def cal_val(self) -> None:
        return self.order.index(self.reg) + self.pos * (10 ** -10)
    
    def __str__(self) -> str:
        return "%s:%s" % (self.reg, self.pos)
    
    def __repr__(self) -> str:
        return "Location(reg=%r, pos=%r)" % (self.reg, self.pos)

    def __lt__(self, other):
        return self.val < other.val
    
    def __le__(self, other):
        return self.val <= other.val
    
    def __eq__(self, other):
        return self.val == other.val
    
    def __ne__(self, other):
        return self.val != other.val
    
    def __gt__(self, other):
        return self.val > other.val
    
    def __ge__(self, other):
        return self.val >= other.val
    
    def __add__(self, pos):
        return Location(
            self.reg,
            self.pos + pos,
            self.order,
            self.name
        )
    
    def __sub__(self, pos):
        return Location(
            self.reg,
            self.pos - pos,
            self.order,
            self.name
        )