class Node: ## node class
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
        self.branchType='node'
        self.length=0.0 ## branch length, recovered from string
        self.height=None ## height, set by traversing the tree, which adds up branch lengths along the way
        self.absoluteTime=None ## branch end point in absolute time, once calibrations are done
        self.parent=None ## reference to parent node of the node
        self.children=[] ## a list of descendent branches of this node
        self.traits={} ## dictionary that will contain annotations from the tree string, e.g. {'posterior':1.0}
        self.index=None ## index of the character designating this object in the tree string, it's a unique identifier for every object in the tree
        self.childHeight=None ## the youngest descendant tip of this node
        self.x=None ## X and Y coordinates of this node, once drawTree() is called
        self.y=None
        ## contains references to all tips of this node
        self.leaves=set() ## is a set of tips that are descended from it

    def is_leaflike(self):
        return False

    def is_leaf(self):
        return False

    def is_node(self):
        return True
