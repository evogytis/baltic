import logging
from .branchLike import BranchLike

logger = logging.getLogger("baltic.Node")


class Node(BranchLike): ## node class

    def __init__(self):
        super().__init__() # run branchLike initializer
        self.branchType='node'
        self.children=[] ## a list of descendent branches of this node
        self.childHeight=None ## the youngest descendant tip of this node
        self.leaves=set() ## is a set of tips that are descended from it

    def is_leaflike(self):
        return False

    def is_leaf(self):
        return False

    def is_node(self):
        return True
