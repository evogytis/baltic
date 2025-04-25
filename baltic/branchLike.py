import logging

logger = logging.getLogger("baltic.BranchLike")


class BranchLike:

    def __init__(self):
        self.branchType=None
        self.length=0.0 ## branch length, recovered from string
        self.height=None ## height, set by traversing the tree, which adds up branch lengths along the way
        self.absoluteTime=None ## branch end point in absolute time, once calibrations are done
        self.parent=None ## reference to parent node of the node
        self.traits={} ## dictionary that will contain annotations from the tree string, e.g. {'posterior':1.0}
        self.index=None ## index of the character designating this object in the tree string, it's a unique identifier for every object in the tree
        self.x=None ## X and Y coordinates of this node, once drawTree() is called
        self.y=None
        ## contains references to all tips of this node

    def get_path_to_root(self, path=None):
        if path is None: path = []

        #TODO: make this self.parent.parent or flag the superroot
        if self.parent is None: # the root doesn't have a parent, so just return the path including the root
            return path + [self]
        else: # if we aren't at the root, recurse
            newPath = path + [self]
            return self.parent.get_path_to_root(newPath)

    def get_siblings(self, include_self=False):
        if self.parent is None:
            logger.warning("Attempted to find siblings of root node. Returning empty list.")
            return []

        sibs = self.parent.children

        if include_self:
            return sibs

        sib_set = set(sibs)
        sib_set.remove(self)
        sibs = list(sib_set)
        return sibs

    def is_leaflike(self):
        return False

    def is_leaf(self):
        return False

    def is_node(self):
        return False
