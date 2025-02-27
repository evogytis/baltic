class Leaf: ## leaf class
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
        self.branchType='leaf'
        self.name=None ## name of tip after translation, since BEAST trees will generally have numbers for taxa but will provide a map at the beginning of the file
        self.index=None ## index of the character that defines this object, will be a unique ID for each object in the tree
        self.length=None ## branch length
        self.absoluteTime=None ## position of tip in absolute time
        self.height=None ## height of tip
        self.parent=None ## parent
        self.traits={} ## trait dictionary
        self.x=None ## position of tip on x axis if the tip were to be plotted
        self.y=None ## position of tip on y axis if the tip were to be plotted

    def is_leaflike(self):
        return True

    def is_leaf(self):
        return True

    def is_node(self):
        return False
