from .branchLike import BranchLike

class Leaf(BranchLike): ## leaf class
    """
    Represents a leaf in a phylogenetic tree.
    
    Attributes:
    branchType (str): The type of branch, default is 'leaf'.
    name (str): The name of the tip after translation, assigned in `make_tree()`.
    index (int): The index of the character that defines this object in the tree string, assigned in `make_tree()`.
    length (float): The length of the branch, assigned in `make_tree()`.
    absoluteTime (float or None): The position of the tip in absolute time, assigned in `setAbsoluteTime()`.
    height (float): The height of the tip, assigned in `traverse_tree()`.
    parent (node): The parent node, assigned in `make_tree()`.
    traits (dict): Dictionary containing traits associated with the leaf, assigned in `make_tree()`.
    x (float or None): The x-coordinate for plotting, default is None.
    y (float or None): The y-coordinate for plotting, default is None.

    Docstring generated with ChatGPT 4o.
    """
    def __init__(self):
        super().__init__() # all of leaf's traits come from BranchLike

    def is_leaflike(self):
        return True

    def is_leaf(self):
        return True

    def is_node(self):
        return False
