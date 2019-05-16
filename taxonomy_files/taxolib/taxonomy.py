"""
Provides classes that represent complete taxonomies, built using components from
the taxacomponents module.
"""


from taxacomponents import Citation, RankTable, Taxon
from taxonvisitor import TaxonVisitor
from taxonvisitors_concrete import PrintTaxonVisitor, CSVTaxonVisitor
from nameresolve import CoLNamesResolver


class TaxonomyError(Exception):
    """
    A basic exception class for reporting errors encountered while working with taxonomies.
    """
    def __init__(self, msg):
        msg = 'Taxonomy error:\n  ' + msg
        Exception.__init__(self, msg)


class TaxonomyBase:
    # Define the "nil" UUID constant as returned by the uuid-osp Postgres module
    # function uuid_nil().
    #NIL_UUID = '00000000-0000-0000-0000-000000000000'
    NIL_UUID = 0

    def __init__(self, taxonomy_id, name='', ismaster=False, citation=None, roottaxon=None):
        self.taxonomy_id = taxonomy_id
        self.name = name
        self.ismaster = ismaster
        self.citation = citation
        self.roottaxon = roottaxon

    def loadFromDB(self, pgcur, taxanum=-1, maxdepth=-1):
        """
        Attempts to load the taxonomy from a taxonomy database, including the full tree
        of taxa.  If taxanum > 0, then only taxanum taxa will be loaded.  If maxdepth > -1,
        the taxa tree will only be traversed to a depth of maxdepth.
        """
        query = """SELECT name, citation_id, ismaster, root_tc_id
            FROM taxonomies
            WHERE taxonomy_id=?"""
        pgcur.execute(query, (self.taxonomy_id,))
        res = pgcur.fetchone()

        if res == None:
            raise TaxonomyError('Taxonomy ID ' + str(self.taxonomy_id) + ' was not found in the database.')

        self.name = res[0]
        self.ismaster = res[2]
        roottc_id = res[3]

        # Create the Citation object.
        self.citation = Citation()
        self.citation.loadFromDB(pgcur, res[1])

        # Get the rank ID and taxonomy ID of the root taxon concept.
        query = """SELECT tc.rank_id, tc.taxonomy_id
            FROM taxon_concepts tc, ranks r
            WHERE tc.tc_id=? AND tc.rank_id=r.rank_id"""
        pgcur.execute(query, (roottc_id,))
        res = pgcur.fetchone()
        rankid = res[0]
        root_taxonomy_id = res[1]

        # Initialize the rank lookup table.
        rankt = RankTable()
        rankt.loadFromDB(pgcur)

        # Load the taxa tree.
        self.roottaxon = Taxon(self.taxonomy_id, rankid, rankt, roottaxo_id = root_taxonomy_id, isroot=True)
        self.roottaxon.loadFromDB(pgcur, roottc_id, taxanum, maxdepth)

    def persist(self):
        """
        Persist the Taxonomy to the database.  This method should be implemented by
        concrete subclasses.
        """
        pass

    def __str__(self):
        tstr = 'name: ' + self.name + '\nID: ' + str(self.taxonomy_id) + '\nmaster: '
        if self.ismaster:
            tstr += 'yes'
        else:
            tstr += 'no'

        return tstr

    def printTaxonomyInfo(self):
        """
        Prints the metadata that describes this taxonomy.
        """
        print '** Taxonomy information **'
        print str(self)
        print str(self.citation)

    def printCSVTaxaTree(self, numtaxa=-1, maxdepth=-1):
        """
        Prints the tree of taxa for this taxonomy in "flat" format as CSV outut.  If
        numtaxa > 0, only the first numtaxa taxa will be printed.  If maxdepth > -1,
        the taxa tree will only be traversed to a depth of maxdepth.
        """
        if numtaxa > 0:
            print '(Only printing first', numtaxa, 'taxa.)'
        if maxdepth > -1:
            print '(Only traversing taxa tree to a depth of ' + str(maxdepth) + '.)'

        csvvisitor = CSVTaxonVisitor(numtaxa, maxdepth)
        csvvisitor.visit(self.roottaxon)

    def printTaxaTree(self, numtaxa=-1, maxdepth=-1):
        """
        Prints the tree of taxa for this taxonomy.  If numtaxa > 0, only the first numtaxa
        taxa will be printed.  If maxdepth > -1, the taxa tree will only be traversed to a
        depth of maxdepth.
        """
        print '** Taxa tree **'
        if numtaxa > 0:
            print '(Only printing first', numtaxa, 'taxa.)'
        if maxdepth > -1:
            print '(Only traversing taxa tree to a depth of ' + str(maxdepth) + '.)'

        ptvisitor = PrintTaxonVisitor(numtaxa, maxdepth)
        ptvisitor.visit(self.roottaxon)

    def printAll(self, numtaxa=-1, maxdepth=-1):
        """
        Prints a text representation of this taxonomy, including the tree of taxa.
        If numtaxa > 0, only the first numtaxa taxa will be printed.  If maxdepth > -1,
        the taxa tree will only be traversed to a depth of maxdepth.
        """
        self.printTaxonomyInfo()
        print
        self.printTaxaTree(numtaxa, maxdepth)


class Taxonomy(TaxonomyBase):
    """
    A class that represents a single taxonomy in the MOL taxonomy database.  Provides methods
    to load a taxonomy from the database and persist a taxonomy to the database.  Can also link
    a taxonomy to the backbone taxonomy.
    """
    def __init__(self, taxonomy_id, name='', ismaster=False, citation=None, roottaxon=None):
        TaxonomyBase.__init__(self, taxonomy_id, name, ismaster, citation, roottaxon)

        # A reference for the backbone taxonomy, which encompasses all other taxonomies.
        # This reference is used if this taxonomy is linked to the backbone taxonomy.
        self.bb_taxonomy = None

    def linkToBackbone(self, pgcur, adjustdepth=True):
        """
        Tries to connect this taxonomy to the backbone taxonomy, creating new nodes
        in the backbone taxonomy, if needed, to link the two together. If adjustdepth
        is True, the depth property of all nodes in the taxonomy are set to match the
        correct depth relative to the root of the backbone taxonomy.  Returns True if
        the linking operation succeeded, False otherwise.
        """
        bb_taxonomy = BackboneTaxonomy(pgcur)
        if bb_taxonomy.linkTaxonomy(self):
            self.bb_taxonomy = bb_taxonomy
            if adjustdepth:
                self.bb_taxonomy.setNodeDepths()
            return True
        else:
            self.bb_taxonomy = None
            return False

    def getBackboneTaxonomy(self):
        """
        Returns a reference to the backbone taxonomy object that links this taxonomy
        to the MOL backbone taxonomy.
        """
        return self.bb_taxonomy

    def persist(self, pgcur, printprogress=False):
        """
        Writes the taxonomy information to the database, if it does not already
        exist.  This includes calling the persist() methods on the Citation and
        Taxon tree associated with this Taxonomy object.
        """
        # First, check if this taxonomy already exists in the database.
        query = """SELECT taxonomy_id
            FROM taxonomies
            WHERE taxonomy_id=? AND ismaster=?"""
        pgcur.execute(query, (self.taxonomy_id, self.ismaster))
        res = pgcur.fetchone()

        if res == None:
            # Write the citation information to the database, if needed.
            citation_id = self.citation.persist(pgcur)
        
            # Create the initial database entry for the taxonomy metadata so that the
            # foreign key constraint for the child taxon concepts can be satisfied.
            query = """INSERT INTO taxonomies
                (taxonomy_id, name, citation_id, ismaster, root_tc_id)
                VALUES (?, ?, ?, ?, ?)"""
            pgcur.execute(query, (self.taxonomy_id, self.name, citation_id, self.ismaster, None))

            # Make sure all taxon concepts, including those from the backbone taxonomy,
            # are persisted to the database.  Use the "nil" UUID as the parent_id for
            # the root of the taxonomy if there is not an existing root entry.
            if self.bb_taxonomy != None:
                self.bb_taxonomy.roottaxon.persist(pgcur, self.NIL_UUID, printprogress,
                        self.roottaxon.depth)
            else:
                self.roottaxon.persist(pgcur, self.NIL_UUID, printprogress, self.roottaxon.depth)

            # Get the ID of the root taxon.
            root_tcid = self.roottaxon.existsInDB(pgcur)

            # Update the taxonomy metadata entry with the root taxon concept ID.
            query = """UPDATE taxonomies
                SET root_tc_id=?
                WHERE taxonomy_id=?"""
            pgcur.execute(query, (root_tcid, self.taxonomy_id))
            pgcur.connection.commit()
        elif printprogress:
            print ('The metadata for taxonomy "' + self.name + '" (ID ' + str(self.taxonomy_id) +
                    ') already exist in the database; no changes were made.')

    def printAll(self, numtaxa=-1, maxdepth=-1):
        """
        Prints a text representation of this taxonomy, including the tree of taxa.
        If numtaxa > 0, only the first numtaxa taxa will be printed.  If maxdepth > -1,
        the taxa tree will only be traversed to a depth of maxdepth.  Unlike the method
        in the base class, this method accounts for the possibility of this taxonomy
        being linked to the backbone taxonomy.
        """
        self.printTaxonomyInfo()
        print
        if self.bb_taxonomy != None:
            self.bb_taxonomy.printTaxaTree(numtaxa, maxdepth)
        else:
            self.printTaxaTree(numtaxa, maxdepth)


class DepthAdjustVisitor(TaxonVisitor):
    """
    Sets the "depth" values for all Taxon objects in a taxa tree, using an initial
    starting depth value.
    """
    def __init__(self, startdepth):
        """
        Assigns startdepth as the "depth" value for the top-level Taxon object.  All
        other "depth" values are calculated relative to startdepth.
        """
        TaxonVisitor.__init__(self)

        self.startdepth = startdepth

    def processTaxon(self, taxon, depth):
        taxon.depth = self.startdepth + depth


class BackboneTaxonomy(TaxonomyBase):
    """
    A special case of Taxonomy that represents the MOL backbone taxonomy.  Provides
    methods to link other taxonomies to the backbone taxonomy.  Does not provide a
    persist() method because the backbone taxonomy metadata are set when the database
    tables are created.
    """
    def __init__(self, pgcur):
        """
        Initialize the backbone Taxonomy object and automatically load it from the
        database, but load only the root node by default.
        """
        self.pgcur = pgcur

        # The ID of the backbone taxonomy is always 1.
        TaxonomyBase.__init__(self, 1)
        self.loadFromDB(pgcur)

    def loadFromDB(self, pgcur, taxanum=-1, maxdepth=0):
        """
        Exactly the same as loadFromDB() from the superclass, except loads only the root
        taxonomy node (i.e., Eukaryota) by default.
        """
        TaxonomyBase.loadFromDB(self, pgcur, taxanum, maxdepth)

    def linkTaxonomy(self, taxonomy):
        """
        Given a Taxonomy object, this method searches for the root taxon concept in the
        database, verifies whether it is already connected to the MOL backbone taxonomy,
        and if not, creates the Taxon objects needed to link it to the backbone taxonomy.
        To do this, the method first checks whether the taxon concept exists in the
        database, and if so, whether it has a parent.  If it has a parent, then it is
        assumed to already be connected to the backbone taxonomy.  Provided that
        taxonomies are correctly added to the database, this assumption should be valid.
        If it is not already linked to the backbone taxonomy, then Catalog of Life is
        used to try to infer the missing taxon nodes that connect the target taxonomy
        to the backbone taxonomy.  If needed, the depth properties of nodes in the target
        taxonomy are adjusted so the values are relative to the depth of the target
        taxonomy's root node in the backbone taxonomy.  If the linking is succesful,
        the method returns True; otherwise, False is returned.
        """
        # Load any parent links to the target taxonomy from the database.
        topnode = self.getLinksFromDB(taxonomy)

        # See if we made it back to the root of the backbone taxonomy.
  
        if topnode.equals(self.roottaxon):
            # We did, so simply link the child of the returned node to our root taxon.
            self.roottaxon.addChild(topnode.children[0])
            success = True
        else:
            # Otherwise, try to use Catalog of Life to fill in any missing links.
            success = self._buildCoLLinks(topnode)

        return success

    def _buildCoLLinks(self, taxon):
        """
        Uses Catalog of Life to fill in missing taxa needed to link the target taxon to the
        MOL backbone taxonomy.  If linking was successful, the target taxon will be connected
        to the backbone root taxon by one or more linking taxa.  Returns True on success;
        False otherwise.
        """
        # Use the Catalog of Life names resolver to try to get higher taxonomy information
        # for the taxon.
        resolver = CoLNamesResolver()
        searchres = resolver.searchCoLForTaxon(taxon, taxon.name.namestr, True)
        if searchres == None:
            return False

        res, sname, srank, authorinfo = searchres

        # Process each parent taxon in the CoL classification, creating a chain of Taxon
        # objects to capture the higher taxonomy.  Because the name resolver search method
        # verifies that the kingdom is correct, we already know that we are connecting the
        # taxonomy to the correct kingdom.
        taxaxml = res.find('./classification')
        # It is important that we use the rank system from the taxonomy (not the backbone)
        # to ensure that rank name lookups retrieve the correct ID.
        tranksys = taxon.ranksys
        ranktable = taxon.rankt
        curnode = self.roottaxon
        for taxonxml in taxaxml:
            namestr = taxonxml.find('name').text
            rankstr = taxonxml.find('rank').text
            child = curnode.createChild(ranktable.getID(rankstr, tranksys), namestr)
            #print child
            curnode = child

        # Link the root of the target taxonomy to the backbone taxonomy.
        curnode.addChild(taxon)

        return True

    def getLinksFromDB(self, taxonomy):
        """
        Determines whether a taxonomy's root taxon is already linked to the MOL backbone
        taxonomy by attempting to follow parent_id links back to the backbone root.
        Returns the top-most node that could be reached by following the links upward.
        """
        # See if the root taxon_concept already has a parent.
        curnode = taxonomy.roottaxon
        parent_id = taxonomy.roottaxon.getParentIDFromDB(self.pgcur)

        # Follow parent links upwards until we reach the root or any other node that
        # has no parent or does not yet exist in the database.
        while parent_id != None and parent_id != self.NIL_UUID:
            # Create the parent node and load it from the database.
            parent = Taxon(curnode.taxonomy_id, curnode.rank_id, curnode.rankt)
            parent.loadFromDB(self.pgcur, parent_id, maxdepth=0)

            parent.addChild(curnode)

            curnode = parent
            parent_id = curnode.getParentIDFromDB(self.pgcur)

        return curnode
    
    def setNodeDepths(self):
        """
        After linking a new taxonomy to the backbone taxonomy, the values of the depth
        properties on the Taxon objects in the target taxonomy are likely to be incorrect.
        This method will visit all nodes and set the correct value of the depth property
        for each node.
        """
        depthvisitor = DepthAdjustVisitor(0)
        depthvisitor.visit(self.roottaxon)

