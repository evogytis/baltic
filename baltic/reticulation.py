class Reticulation: ## reticulation class (recombination, conversion, reassortment)
    """
    Represents a reticulation event in a phylogenetic tree, such as recombination or reassortment.
    
    Attributes:
    branchType (str): The type of branch, default is 'leaf'.
    length (float): The length of the branch, assigned in `make_tree()`.
    height (float): The height of the branch, assigned in `traverse_tree()`.
    absoluteTime (float or None): The absolute time of the event, default is None, typically assigned in `setAbsoluteTime()`.
    parent (node): The parent node, assigned in `make_tree()`.
    traits (dict): Dictionary of traits associated with the reticulation event, default is empty.
    index (int): The index of the node in the tree string, assigned in `make_tree()`.
    name (str): The name of the reticulation event.
    x (float or None): The x-coordinate for plotting, default is None.
    y (float or None): The y-coordinate for plotting, default is None.
    width (float): The width of the node for plotting, default is 0.5.
    target (node or leaf): The target node where reticulation event lands, assigned in `make_tree()`.
    
    Docstring generated with ChatGPT 4o.
    """
    def __init__(self,name):
        self.branchType='leaf'
        self.length=0.0
        self.height=0.0
        self.absoluteTime=None
        self.parent=None
        self.traits={}
        self.index=None
        self.name=name
        self.x=None
        self.y=None
        self.width=0.5
        self.target=None

    def is_leaflike(self):
        return True

    def is_leaf(self):
        return False

    def is_node(self):
        return False
