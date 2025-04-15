import logging
from .branchLike import BranchLike

logger = logging.getLogger("baltic.Reticulation")


class Reticulation(BranchLike): ## reticulation class (recombination, conversion, reassortment)

    def __init__(self,name):
        super().__init__() # Inherit traits from BranchLike
        self.branchType='leaf'
        self.name=name
        self.width=0.5
        self.target=None

    def is_leaflike(self):
        return True

    def is_leaf(self):
        return False

    def is_node(self):
        return False
