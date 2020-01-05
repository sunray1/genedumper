"""
Concrete child class implementations of TaxonVisitor.  All of these classes
perform operations on a taxonomic tree by recursively traversing the tree and
visiting all taxa within the tree.
This module is not very well named, nor is the code here particularly well
organized.  The classes here should probably be organized into modules by
function, similar to the way the TaxonVisitor child classes for name resolution
are organized.
"""

from taxolib.taxonvisitor import TaxonVisitor
import sys
from taxolib.csvtaxonomy import UnicodeDictWriter


class PrintTaxonVisitor(TaxonVisitor):
    """
    Prints a simple text representation of a taxon tree.
    """
    def processTaxon(self, taxon, depth):
        print('   ' * depth + taxon.__str__())
        synstr = taxon.getSynonymsString()
        if synstr != '':
            print('   ' * depth + '  -- Synonyms: ' + synstr)


class RankAccumulatorTaxonVisitor(TaxonVisitor):
    """
    A taxon visitor class that builds a list of all taxonomic ranks that are used in
    a taxonomy tree.
    """
    def visit(self, taxon):
        """
        Returns a list of all taxonomic rank names used in a taxonomy tree.  The rank
        names will be sorted according to their IDs in the database.
        """
        # Initialize a list for the rank strings.
        self.ranks = []

        # Call the superclass method implementation to initiate the traversal.
        TaxonVisitor.visit(self, taxon)

        # Sort the ranks according to their rank IDs.
        self.ranks.sort(key=lambda rankstr: taxon.rankt.getID(rankstr, taxon.ranksys))

        return self.ranks

    def processTaxon(self, taxon, depth):
        rankstr = taxon.getRankString()

        if rankstr not in self.ranks:
            self.ranks.append(rankstr)


class CSVTaxonVisitor(TaxonVisitor):
    """
    Prints a CSV representation of a taxon tree.
    """
    def visit(self, taxon):
        # Get a list of all ranks used in the tree..
        ranks = RankAccumulatorTaxonVisitor().visit(taxon)

        # Set up the CSV writer.
        csvheader = ranks + ['Author', 'Synonyms', 'Citation']
        self.writer = UnicodeDictWriter(sys.stdout, fieldnames=csvheader)
        self.writer.writeheader()

        # Create a dictionary for storing CSV row values.
        self.rowvals = {}
        for colname in csvheader:
            self.rowvals[colname] = ''

        # Call the superclass method implementation to initiate the traversal.
        TaxonVisitor.visit(self, taxon)

        return

    def processTaxon(self, taxon, depth):
        rank = taxon.getRankString()
        self.rowvals[rank] = taxon.name.namestr

        # Only get author and synonym information and print a row if we're at a
        # leaf of the tree.  This makes the CSV output consistent with the way
        # CSV taxonomy files are parsed by the input routines.
        if len(taxon.children) == 0:
            self.rowvals['Author'] = taxon.getAuthorDisplayString()
            self.rowvals['Synonyms'] = taxon.getSynonymsString()
            self.rowvals['Citation'] = taxon.name.getCitationString()
            self.writer.writerow(self.rowvals)

    def postTaxonProcessing(self, taxon, depth):
        self.rowvals[taxon.getRankString()] = ''
        self.rowvals['Author'] = ''
        self.rowvals['Synonyms'] = ''
        self.rowvals['Citation'] = ''


class NamesTaxonVisitor(TaxonVisitor):
    """
    A taxon visitor class that builds a list of all Name objects in the taxa tree.
    The results are returned as a list of tuples, where each tuple contains the
    Name object and a string representing the taxonomic rank of the name; i.e.:
    (Name, rankstring).
    """
    def visit(self, taxon):
        # Initialize a list for the Name objects.
        self.names = []

        # Call the superclass method implementation.
        TaxonVisitor.visit(self, taxon)

        return self.names

    def processTaxon(self, taxon, depth):
        self.names.append((taxon.name, taxon.getRankString()))

