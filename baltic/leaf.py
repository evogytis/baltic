import logging
from .branchLike import BranchLike

logger = logging.getLogger("baltic.Leaf")


class Leaf(BranchLike): ## leaf class

    def __init__(self,name):
        super().__init__() # all of leaf's traits come from BranchLike
        self.branchType='leaf'
        self.name=name

    def is_leaflike(self):
        return True

    def is_leaf(self):
        return True

    def is_node(self):
        return False
