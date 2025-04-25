import logging
from .branchLike import BranchLike

logger = logging.getLogger("baltic.Clade")


class Clade(BranchLike): ## clade class

    def __init__(self,givenName):
        self.branchType='leaf' ## clade class poses as a leaf
        self.subtree=None ## subtree will contain all the branches that were collapsed
        self.leaves=None
        self.name=givenName ## the pretend tip name for the clade
        self.lastHeight=None ## refers to the height of the highest tip in the collapsed clade
        self.lastAbsoluteTime=None ## refers to the absolute time of the highest tip in the collapsed clade
        self.width=1

    def is_leaflike(self):
        return True

    def is_leaf(self):
        return False

    def is_node(self):
        return False
