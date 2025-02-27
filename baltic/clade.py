class Clade: ## clade class
    """
    Represents a collapsed clade in a phylogenetic tree.
    
    Attributes:
    branchType (str): The type of branch, default is 'leaf'.
    subtree (list): The subtree containing all the branches that were collapsed, assigned in `collapseSubtree()`.
    leaves (set): Names of descendant tips in the collapsed clade, assigned in `collapseSubtree()`.
    length (float): The length of the branch, assigned in `collapseSubtree()`.
    height (float or None): The height of the branch, assigned in `collapseSubtree()`.
    absoluteTime (float or None): The absolute time of the event, assigned in `collapseSubtree()`.
    parent (node): The parent node, assigned in `collapseSubtree()`.
    traits (dict): Dictionary of traits associated with the clade, assigned in `collapseSubtree()`.
    index (int or None): The index of the node, assigned in `collapseSubtree()`.
    name (str): The name assigned to the clade when it was collapsed.
    x (float or None): The x-coordinate for plotting, default is None.
    y (float or None): The y-coordinate for plotting, default is None.
    lastHeight (float): The height of the highest tip in the collapsed clade, assigned in `collapseSubtree()`.
    lastAbsoluteTime (float or None): The absolute time of the highest tip in the collapsed clade, assigned in `collapseSubtree()`.
    width (float): The width of the node for plotting, default is 1.
    
    Docstring generated with ChatGPT 4o.
    """
    def __init__(self,givenName):
        self.branchType='leaf' ## clade class poses as a leaf
        self.subtree=None ## subtree will contain all the branches that were collapsed
        self.leaves=None
        self.length=0.0
        self.height=None
        self.absoluteTime=None
        self.parent=None
        self.traits={}
        self.index=None
        self.name=givenName ## the pretend tip name for the clade
        self.x=None
        self.y=None
        self.lastHeight=None ## refers to the height of the highest tip in the collapsed clade
        self.lastAbsoluteTime=None ## refers to the absolute time of the highest tip in the collapsed clade
        self.width=1

    def is_leaflike(self):
        return True

    def is_leaf(self):
        return False

    def is_node(self):
        return False
