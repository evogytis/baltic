from .branchLike import BranchLike

class Node(BranchLike): ## node class
    """
    Represents a node in a phylogenetic tree.
    
    Attributes:
    branchType (str): The type of branch, default is 'node'.
    length (float): The length of the branch, assigned in `make_tree()`.
    height (float): The height of the branch, assigned in `traverse_tree()`.
    absoluteTime (float or None): The branch endpoint in absolute time, assigned in `setAbsoluteTime()`.
    parent (node): The parent node, assigned in `make_tree()`.
    children (list): A list of descendant branches of this node, assigned in `make_tree()`.
    traits (dict): Dictionary containing annotations from the tree string, assigned in `make_tree()`.
    index (int): The index of the character designating this object in the tree string, a unique identifier.
    childHeight (float or None): The height of the youngest (last) descendant tip of this node, assigned in `traverse_tree()`.
    x (float or None): The x-coordinate for plotting, default is None.
    y (float or None): The y-coordinate for plotting, default is None.
    leaves (set): A set of tips that are descended from this node, assigned in `traverse_tree()`.
    
    Docstring generated with ChatGPT 4o.
    """


    def __init__(self):
        super().__init__() # run branchLike initializer
        self.branchType='node'
        self.childHeight=None ## the youngest descendant tip of this node
        self.leaves=set() ## is a set of tips that are descended from it

    def is_leaflike(self):
        return False

    def is_leaf(self):
        return False

    def is_node(self):
        return True
