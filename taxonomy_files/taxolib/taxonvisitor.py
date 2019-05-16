"""
Provides a base class that implements operations on a taxonomic tree that
require recursively visiting each taxon within the tree.
"""


class TaxonVisitor:
    """
    Base class for all taxon tree visitor classes.  Client code calls the visit() method
    with an instance of Taxon as an argument.  The visitor class will then traverse the
    taxon tree and operate on each Taxon object in the tree.  The traversal can be limited
    by either the total number of taxa processed or tree depth (or both).  In essence, this
    base class encapsulates an algorithm for traversing a taxa tree and allows operations
    on tree objects to be implemented independently of the Taxon implementation.
    """
    def __init__(self, numtaxa=-1, maxdepth=-1):
        """
        If numtaxa > 0, only the first numtaxa taxa will be visited when a taxon tree
        is traversed.  If maxdepth > -1, the tree will only be traversed to a depth
        of maxdepth.
        """
        self.numtaxa = numtaxa
        self.maxdepth = maxdepth

    def visit(self, taxon):
        """
        Initiates the taxon tree traversal.
        """
        self.taxacnt = 0

        self._traverseTree(taxon, 0)

    def _traverseTree(self, taxon, depth):
        """
        Internal method for traversing a taxon tree that tracks the recursion depth.
        """
        self.processTaxon(taxon, depth)

        self.taxacnt += 1

        if (self.maxdepth < 0 or depth < self.maxdepth) and self.doRecursion():
            for child in taxon.children:
                if self.numtaxa > 0 and self.taxacnt >= self.numtaxa:
                    break
                self._traverseTree(child, depth + 1)

        self.postTaxonProcessing(taxon, depth)

    def doRecursion(self):
        """
        This method can be overriden by child classes and used to implement additional
        criteria for deciding whether to recursively descend into the next level of a
        taxa tree.  The method is called prior to processing the child taxa of a taxon;
        if it returns True, the recursion is continued, otherwise, the children are
        not visited.
        """
        return True

    def processTaxon(self, taxon, depth):
        """
        This method is called for each Taxon object in the tree.  The argument 'depth'
        provides the current depth in the tree, with the root at 0.  This method should
        be overridden by child classes to actually do something with each Taxon object.
        """
        pass

    def postTaxonProcessing(self, taxon, depth):
        """
        This method is called after tree traversal has returned from recursively
        traversing taxon's descendents.  It can be overridden by child classes to
        implement "clean up" code that should be run before leaving a taxon.
        """
        pass

