"""
Provides classes that model the basic pieces of a taxonomy: rank information, taxon concepts,
names, citations, and algorithms to operate on trees of taxon concepts.
"""


class RankTable:
    """
    Provides simple rank name or rank ID lookup (i.e., lookup one, return the other).  Rank
    ID lookups are always unambiguous, but rank name lookups require a rank system ID.
    RankTable keeps track of a default rank system to use for name lookups.
    """
    def __init__(self):
        # Define two dictionaries so we can efficiently search in both directions.
        self.IDstonames = {}
        self.namestoIDs = {}

    def loadFromDB(self, pgcur):
        """
        Loads the rank ID and name information from the MOL taxonomy database.
        """
        # Get all of the rank system IDs.
        pgcur.execute('SELECT ranksys_id FROM rank_systems')
        for rec in pgcur:
            self.namestoIDs[rec[0]] = {}

        # Now get the ranks for each rank system.
        for rsid in self.namestoIDs.keys():
            pgcur.execute('SELECT rank_id, namestr FROM ranks WHERE ranksys_id=? ORDER BY rank_id',
                    (rsid,))
            for rec in pgcur:
                self.namestoIDs[rsid][rec[1]] = rec[0]
                self.IDstonames[rec[0]] = rec[1]

    def getMinRankID(self, ranksys):
        """
        Returns the minimum ID within a given rank system.
        """
        return min(self.namestoIDs[ranksys].itervalues())

    def getMaxRankID(self, ranksys):
        """
        Returns the minimum ID within a given rank system.
        """
        return max(self.namestoIDs[ranksys].itervalues())

    def isIDInRanksys(self, ranksys, rankid):
        """
        Checks if a given rank ID is within a given rank system.
        """
        return rankid in self.namestoIDs[ranksys].itervalues()

    def getRanksysByRankID(self, rankid):
        """
        Gets the rank system to which the given rankid belongs.  Returns None if the
        lookup failed.
        """
        for ranksys in self.namestoIDs.keys():
            try:
               if rankid in self.namestoIDs[ranksys].itervalues():
                   return ranksys
            except:
               if rankid in self.namestoIDs[ranksys].items():
                   return ranksys

        return None

    def getName(self, rankid):
        """
        Lookup a rank ID and return the associated name string.
        """
        return self.IDstonames[rankid]

    def getID(self, namestr, ranksys):
        """
        Lookup a rank name string for a given rank system and return the associated
        rank ID.
        """
        return self.namestoIDs[ranksys][namestr]


class NameList:
    """
    Represents a set of names that are associated with a particular taxon concept.
    """
    def __init__(self):
        # The list of Name objects.
        self.names = []
        # A corresponding list indicating whether to use parentheses for each name author string.
        self.useparenslist = []

        # The index of the preferred (valid) name in the names list.
        self.preferred = None

        # The taxon associated with the list of names.
        self.taxon = None

    def getPreferred(self):
        if self.preferred != None:
            return self.names[self.preferred]
        else:
            return None

    def getSynNameStrs(self):
        """
        Returns a list of all name strings that are considered synonyms (i.e., not preferred).
        """
        namelist = []

        for cnt in range(len(self.names)):
            if cnt != self.preferred:
                namelist.append(self.names[cnt].namestr)

        return namelist

    def associateWith(self, taxon):
        """
        Associates this list of names with a particular Taxon object.
        """
        self.taxon = taxon

    def addName(self, name, preferred, useparens):
        """
        Adds a Name object to the list of names.  The value of preferred indicates whether
        the new name is the preferred name for the taxon.
        """
        self.names.append(name)
        self.useparenslist.append(useparens)

        if preferred:
            self.preferred = len(self.names) - 1

    def getPreferredUseParens(self):
        """
        Whether citations for the preferred name should use parentheses.
        """
        return self.useparenslist[self.preferred]

    def setPreferredUseParens(self, useparens):
        """
        Sets whether citations for the preferred name should use parentheses.
        """
        self.setUseParens(self.preferred, useparens)

    def setUseParens(self, index, useparens):
        """
        Sets whether citations for the name at index should use parentheses.
        """
        self.useparenslist[index] = useparens

    def loadFromDB(self, pgcur, tc_id):
        """
        Loads the name list for a given taxon ID from the database.
        """
        self.names = []
        self.useparenslist = []

        query = """SELECT nttc.name_id, nttc.validity, nttc.authordisp_prefix
            FROM names_to_taxonconcepts nttc
            WHERE nttc.tc_id=?"""
        pgcur.execute(query, (tc_id,))
        results = pgcur.fetchall()

        for res in results:
            nameid = res[0]
            if nameid != None:
                # Load the name (and citation information).
                name = Name()
                name.loadFromDB(pgcur, nameid)
                preferred = (res[1] == 'valid')
                use_parens = (res[2] == '(')

                self.addName(name, preferred, use_parens)

    def persist(self, pgcur, tc_id, do_commit=True):
        """
        Ensures that all names in the list are associated with the taxon concept in the
        taxonomy database.  Note that all names besides the preferred/valid name are
        entered into the database as synonyms.
        """
        for name, useparens, index in zip(self.names, self.useparenslist, range(len(self.names))):
            dbwrite = False

            # Process the name and get the name_id.
            name_id = name.persist(pgcur, do_commit)

            # See if this name is already linked to the taxon concept.
            query = """SELECT name_id
                FROM names_to_taxonconcepts nttc
                WHERE nttc.name_id=? AND nttc.tc_id=?"""
            pgcur.execute(query, (name_id, tc_id))

            # If not, link the name to the new taxon_concept.
            if pgcur.fetchone() == None:
                query = """INSERT INTO names_to_taxonconcepts
                    (tc_id, name_id, validity, authordisp_prefix, authordisp_postfix)
                    VALUES (?, ?, ?, ?, ?)"""
                if useparens:
                    pgcur.execute(query, (tc_id, name_id, 'valid' if (index == self.preferred) else 'synonym', '(', ')'))
                else:
                    pgcur.execute(query, (tc_id, name_id, 'valid' if (index == self.preferred) else 'synonym', '', ''))
                dbwrite = True

        if do_commit and dbwrite:
            # End the transaction.
            pgcur.connection.commit()


class Taxon:
    """
    Represents a single taxon concept that approximately models a taxon concept in the MOL
    taxonomies database.  The key difference is that each Taxon keeps a list of references
    to all of its children (rather than a reference to its parent, as in the database).
    """
    @staticmethod
    def find(pgcur, searchstr, ranktable, taxonomy_id=None, rank_id=None, pref_names_only=False):
        """
        A static method to search for taxa in the taxonomy database that have a name matching
        the value of searchstr.  Matching taxa are returned as a list of Taxon objects.  Taxa
        searches can be further limited by providing a taxonomy ID or a rank ID.  If
        pref_names_only == True, only preferred names are searched.
        """
        taxalist = []

        query = """SELECT DISTINCT tc.tc_id, tc.taxonomy_id, tc.rank_id
            FROM names n, names_to_taxonconcepts nttc, taxon_concepts tc
            WHERE n.namestr LIKE ? AND n.name_id=nttc.name_id AND nttc.tc_id=tc.tc_id"""
        params = [searchstr]
        if pref_names_only:
            query += ' AND nttc.preferred=TRUE'
        if taxonomy_id != None:
            query += ' AND tc.taxonomy_id=?'
            params.append(taxonomy_id)
        if rank_id != None:
            query += ' AND tc.rank_id=?'
            params.append(rank_id)
        pgcur.execute(query, params)

        for res in pgcur.fetchall():
            taxon = Taxon(res[1], res[2], ranktable)
            taxon.loadFromDB(pgcur, res[0], maxdepth=0)
            taxalist.append(taxon)

        return taxalist

    def __init__(self, taxonomy_id, rank_id, ranktable, depth=0, namestr='', roottaxo_id=None, isroot=False):
        """
        Most of the method parameters should be self-explanatory.  If this taxon is the root
        of a taxon tree, it can use a different taxonomy ID from its children.  In this case,
        specify isroot = True and provide the root taxonomy ID as the value of roottaxo_id.
        Regardless of whether an alternative taxonomy ID is used for the root, all descendents
        will, by default, use the value of taxonomy_id for their taxonomy ID.
        """
        self.taxonomy_id = taxonomy_id
        self.tc_id = 0
        self.rankt = ranktable
        self.depth = depth
        self.isroot = isroot
        if roottaxo_id == None:
            self.roottaxo_id = self.taxonomy_id
        else:
            self.roottaxo_id = roottaxo_id

        self.namelist = NameList()
        if namestr != '':
            name = Name(namestr)
            self.namelist.addName(name, True, False)

        self.children = []

        self.setRankID(rank_id)

    def __getattr__(self, attrname):
        """
        Implements a dynamic, computed attribute called "name" that points directly to the
        preferred name for this Taxon.  This is merely a convenience for client code.
        """
        if attrname == 'name':
            return self.namelist.getPreferred()
        else:
            raise AttributeError("Taxon instance has no attribute '" + attrname + "'")

    def equals(self, taxon):
        """
        Tests whether two taxon objects are equivalent.  If they have the same taxonomy ID,
        rank ID, depth, and name string, they are considered to be the same.
        """
        return (self.taxonomy_id == taxon.taxonomy_id and self.rank_id == taxon.rank_id and
                self.depth == taxon.depth and self.name.namestr == taxon.name.namestr)

    def setRankID(self, rank_id):
        """
        Sets the rank ID for this Taxon and uses the rank ID to set the correct rank system.
        """
        self.rank_id = rank_id

        # Set which rank system to use.
        #print self.rankt.ranksys, self.rank_id
        #print self.rankt.isIDInRanksys(self.rank_id)
        self.ranksys = self.rankt.getRanksysByRankID(self.rank_id)

    def __str__(self):
        authorstr = self.getAuthorDisplayString()
        if authorstr != '':
            authorstr = ' ' + authorstr

        return self._getRankNameString() + authorstr + '; depth: ' + str(self.depth)

    def getAuthorDisplayString(self):
        """
        Returns the author display string for the valid name for this taxon.
        """
        authorstr = ''
        cite = self.name.getCitation()
        # Get the author display string, adding parentheses if needed.
        if cite != None:
            if cite.authordisp != '' and cite.authordisp != None:
                if self.namelist.getPreferredUseParens():
                    authorstr = '(' + cite.authordisp + ')'
                else:
                    authorstr = cite.authordisp

        return authorstr

    def _getRankNameString(self):
        return self.getRankString() + ': ' + self.name.namestr

    def getRankString(self):
        return self.rankt.getName(self.rank_id)

    def getSynonymsString(self):
        synlist = self.namelist.getSynNameStrs()
        if len(synlist) != 0:
            return ', '.join(synlist)
        else:
            return ''

    def getUseParens(self):
        """
        Reports whether or not parentheses should be used with author display strings for
        citations for the preferred name.
        """
        return self.namelist.getPreferredUseParens()

    def setUseParens(self, useparens):
        """
        Specify whether or not parentheses should be used with author display strings for
        citations for the preferred name.
        """
        self.namelist.setPreferredUseParens(useparens)

    def loadFromDB(self, pgcur, tc_id, taxanum=-1, maxdepth=-1):
        """
        Attempts to load this taxon from the database using the specified tc_id.  All
        children and descendents of this taxon will also be recursively loaded.  Use
        taxanum to limit the number of taxa that are loaded.  If taxanum < 1, all taxa
        will be loaded.  If maxdepth > -1, the tree will only be traversed to a depth
        of maxdepth.
        """
        self._loadFromDB(pgcur, tc_id, taxanum, 0, maxdepth, 0)

    def _loadFromDB(self, pgcur, tc_id, taxanum, taxacnt, maxdepth, curdepth):
        """
        An internal method to recursively load a taxa tree from the taxonomy database.
        Keeps track of how many taxa have been loaded and the depth of the recursion.
        """
        # Get the taxon concept information.
        query = """SELECT tc.rank_id, tc.depth, tc.taxonomy_id
            FROM taxon_concepts tc
            WHERE tc.tc_id=?"""
        pgcur.execute(query, (tc_id,))
        res = pgcur.fetchone()

        self.setRankID(res[0])

        self.tc_id = tc_id
        self.depth = res[1]

        if self.isroot:
            self.roottaxo_id = res[2]
        else:
            self.taxonomy_id = res[2]

        # Load the names (and citation information).
        self.namelist.loadFromDB(pgcur, tc_id)

        taxacnt += 1

        # Now process all children of this taxon that are part of the same taxonomy.
        query = """SELECT tc_id
            FROM taxon_concepts
            WHERE parent_id=? AND taxonomy_id=?"""
        pgcur.execute(query, (tc_id, self.taxonomy_id))
        results = pgcur.fetchall()
        
        if maxdepth < 0 or curdepth < maxdepth:
            for res in results:
                if taxanum > 0 and taxacnt >= taxanum:
                    break
                childtaxon = Taxon(self.taxonomy_id, self.rank_id, self.rankt)
                self.children.append(childtaxon)
                taxacnt = childtaxon._loadFromDB(pgcur, res[0], taxanum, taxacnt, maxdepth, curdepth + 1)

        return taxacnt

    def findChild(self, rank_id, namestr):
        """
        Search for a child taxon by rank ID and preferred name string.  If a match is
        found, the child Taxon object is returned.  Otherwise, None is returned.
        """
        found = False
        index = 0
        while not(found) and index < len(self.children):
            if self.children[index].rank_id == rank_id and self.children[index].name.namestr == namestr:
                found = True
            index += 1
        
        if found:
            return self.children[index - 1]
        else:
            return None

    def addChild(self, taxon):
        """
        Adds a Taxon object as a child of this Taxon object.
        """
        self.children.append(taxon)

    def createChild(self, rank_id, namestr):
        """
        Creates a new Taxon object and adds it as a child of this Taxon object.
        The child object will have the same taxonomy ID and rank lookup table as
        the parent, unless the parent is the root of a tree and has an alternative
        root taxonomy ID.  The child's depth will be self.depth + 1.  The child
        Taxon object is the return value.
        """
        newchild = Taxon(self.taxonomy_id, rank_id, self.rankt, self.depth + 1, namestr)
        self.children.append(newchild)

        return newchild

    def existsInDB(self, pgcur):
        """
        Checks whether this Taxon already exists in the taxonomy database.  This search
        assumes that for a given taxonomy, there are no taxon concepts with the same
        preferred name with the same taxonomic rank.  This should certainly be a reasonable
        assumption (e.g., two families in the same taxonomy cannot have the same name).
        If a matching taxon concept is found, the taxon concept's ID is returned.
        Otherwise, None is returned.
        """
        # Choose which taxonomy ID to use depending on whether this is a root taxon.
        if self.isroot:
            taxo_id = self.roottaxo_id
        else:
            taxo_id = self.taxonomy_id

        # See if a corresponding taxon concept already exists in the database.
        query = """SELECT tc.tc_id
            FROM taxon_concepts tc, names_to_taxonconcepts nttc, names n
            WHERE tc.tc_id=nttc.tc_id AND nttc.name_id=n.name_id AND nttc.validity='valid'
                AND tc.taxonomy_id=? AND tc.rank_id=? AND n.namestr=?"""
        pgcur.execute(query, (taxo_id, self.rank_id, self.name.namestr))
        res = pgcur.fetchone()

        if res == None:
            return None
        else:
            return res[0]

    def getParentIDFromDB(self, pgcur):
        """
        Attempts to retrieve the parent ID for this Taxon from the database.  Returns
        None on failure.
        """
        tc_id = self.existsInDB(pgcur)
        if tc_id == None:
            return  None

        query = """SELECT parent_id
            FROM taxon_concepts tc
            WHERE tc.tc_id=?"""
        pgcur.execute(query, (tc_id,))
        res = pgcur.fetchone()

        return res[0]

    def persist(self, pgcur, parent_id, printprogress=False, rootdepth=0):
        """
        If the taxon_concept corresponding with this Taxon object does not exist in the
        database, this method writes it to the database.  All descendent taxa of this Taxon
        object are then processed recursively.  If printprogress == True, occasional progress
        updates are printed to standard out.  Returns the tc_id of either the existing
        taxon_concept or the newly created taxon_concept.
        This method keeps track of the root depth so it can track its position relative
        to the root node used by the initial persist() call.
        """
        # Choose which taxonomy ID to use depending on whether this is a root taxon.
        if self.isroot:
            taxo_id = self.roottaxo_id
        else:
            taxo_id = self.taxonomy_id

        if printprogress and (self.depth - rootdepth) == 1:
            print('Processing ' + self._getRankNameString() + ' and its descendents...')

        tc_id = self.existsInDB(pgcur)

        if tc_id == None:
            # Add a new taxon concept to the database.
            query = """INSERT INTO taxon_concepts
                (parent_id, taxonomy_id, rank_id, depth)
                VALUES (?, ?, ?, ?)"""
            pgcur.execute(query, (parent_id, taxo_id, self.rank_id, self.depth))
            tc_id = pgcur.lastrowid
        else:
            print ('The taxon concept "' + self._getRankNameString() + '" for taxonomy ID '
                + str(taxo_id) + ' already exists in the database.')

        # Make sure the names are linked to the taxon concept.
        self.namelist.persist(pgcur, tc_id, False)

        # End the transaction.
        pgcur.connection.commit()

        # Process all children of this Taxon object.
        for child in self.children:
            child.persist(pgcur, tc_id, printprogress, rootdepth)

        return tc_id


class Name:
    def __init__(self, namestr='', citation=None):
        self.namestr = namestr
        self.citation = citation
        self.idnum = None

    def setCitation(self, citation):
        self.citation = citation

    def getCitation(self):
        return self.citation

    def getCitationString(self):
        """
        Returns the full citation string for this name, or the empty string
        if no citation information is available.
        """
        if self.citation != None:
            return self.citation.citestr
        else:
            return ''

    def updateCitation(self, citestr, authordisp):
        """
        Augments the existing citation information for a name or creates a new
        Citation object and adds it to the name.  This method currently only
        handles the two text citation fields.  This method will not overrite
        the values of fields that already contain citation data.
        """
        if citestr == None:
            citestr = ''
        if authordisp == None:
            authordisp = ''

        cite = self.getCitation()
        if cite != None:
            # Do not erase existing citation data.
            if (cite.citestr == None or cite.citestr == '') and citestr != '':
                cite.citestr = citestr
            if (cite.authordisp == None or cite.authordisp == '') and authordisp != '':
                cite.authordisp = authordisp
        elif citestr != '' or authordisp != '':
            # Only create a new citation if we have at least some citation data.
            newcite = Citation(citestr, authordisp)
            self.setCitation(newcite)

    def loadFromDB(self, pgcur, name_id):
        """
        Populate this Name object by loading the name data from the database,
        given a unique name_id.
        """
        query = """SELECT name_id, namestr, citation_id
            FROM names
            WHERE name_id=?"""
        pgcur.execute(query, (name_id,))
        res = pgcur.fetchone()

        self.idnum = res[0]
        self.namestr = res[1]
        citeid = res[2]
        if citeid != None:
            self.citation = Citation()
            self.citation.loadFromDB(pgcur, citeid)

    def loadFromDBbyName(self, pgcur, taxonomy_id, rank_id, namestr=''):
        """
        Populate this Name object by loading the name data from the database.
        This method searches for a particular name string, assigned to a particular
        rank, as part of a taxonomy. If no name string is provided, the value of
        this Name's existing name string is used.  The query results are considered
        a match if a single, unique result is returned.  Because this method does
        not use either author information or the name_id, there is always a risk
        that an incorrect name will be loaded (that is, one for which the name string
        is correct but the authorship is incorrect; a homonym of the intended name).
        To help mitigate this, this method will only consider names that are marked
        as "preferred", that are part of a specified target taxonomy, and that are
        associated with a particular taxonomic rank.  This should all make it reliable
        in theory, because the nomenclature codes do not allow multiple taxa with the
        same name at the same rank.
        """
        if namestr == '':
            namestr = self.namestr

        query = """SELECT DISTINCT n.name_id, n.citation_id
            FROM names n, names_to_taxonconcepts nttc, taxon_concepts tc
            WHERE n.namestr=? AND n.name_id=nttc.name_id AND nttc.tc_id=tc.tc_id
                AND nttc.preferred=TRUE AND tc.taxonomy_id=? AND tc.rank_id=?"""
        pgcur.execute(query, (namestr, taxonomy_id, rank_id))
        results = pgcur.fetchall()

        # Check if we only found a single (name_id, citation_id) combination.  If
        # we got multiple matches, we don't know which is correct.
        if results != None and len(results) == 1:
            res = results[0]
            self.namestr = namestr
            self.idnum = res[0]
            citeid = res[1]
            if citeid != None:
                self.citation = Citation()
                self.citation.loadFromDB(pgcur, citeid)

    def persist(self, pgcur, do_commit=True):
        """
        If this name does not exist in the database, creates a new name entry in the
        database.  The name_id of either the existing name or newly created name is returned.
        """
        # First, process the citation and get the citation ID.
        if self.citation != None:
            cite_id = self.citation.persist(pgcur, do_commit)
        else:
            cite_id = None

        # See if there is already an entry for this name in the database.
        # If a name_id has already been retrieved from the database, there is
        # no need for further database querying.
        if self.idnum != None:
            return self.idnum

        query = """SELECT n.name_id
            FROM names n
            WHERE n.namestr=? AND n.citation_id"""
        if cite_id != None:
            pgcur.execute(query + '=?', (self.namestr, cite_id))
        else:
            pgcur.execute(query + ' IS NULL', (self.namestr,))
        res = pgcur.fetchone()

        if res == None:
            # Add the name to the database.
            query = """INSERT INTO names
                (namestr, citation_id)
                VALUES (?, ?)"""
            pgcur.execute(query, (self.namestr, cite_id))

            name_id = pgcur.lastrowid
        else:
            # Use the ID of the existing name entry.
            name_id = res[0]

        if do_commit:
            # End the transaction.
            pgcur.connection.commit()

        return name_id


class Citation:
    def __init__(self, citestr='', authordisp='', doi=None, url=None):
        self.citestr = citestr
        self.authordisp = authordisp
        self.doi = doi
        self.url = url

    def __str__(self):
        cstr = self.citestr
        if self.doi != None:
            cstr += '\nDOI: ' + self.doi
        if self.url != None:
            cstr += '\nURL: ' + self.url

        return cstr
    
    def loadFromDB(self, pgcur, cite_id):
        query = """SELECT citationstr, url, doi, authordisplay
            FROM citations
            WHERE citation_id=?"""
        pgcur.execute(query, (cite_id,))
        res = pgcur.fetchone()

        self.citestr = res[0]
        self.authordisp = res[3]
        self.doi = res[2]
        self.url = res[1]

    def persist(self, pgcur, do_commit=True):
        """
        Writes this citation to the database if it does not already exist.  Returns the
        ID of either the existing citation or the newly created citation.
        """
        # See if matching citation information already exists in the database.
        # Begin by attempting to match the DOI, then the URL, then the full
        # citation string, then the author display string.  Matches based on the
        # author display string only are necessarily "fuzzy", so for an existing
        # record to match, it must also not contain any other citation information.
        citation_id = None
        if self.doi != None and self.doi != '':
            query = """SELECT citation_id
                FROM citations c
                WHERE c.doi=?"""
            pgcur.execute(query, (self.doi,))
            res = pgcur.fetchone()
            if res != None:
                citation_id = res[0]
        elif self.url != None and self.url != '':
            query = """SELECT citation_id
                FROM citations c
                WHERE c.url=?"""
            pgcur.execute(query, (self.url,))
            res = pgcur.fetchone()
            if res != None:
                citation_id = res[0]
        elif self.citestr != None and self.citestr != '':
            query = """SELECT citation_id
                FROM citations c
                WHERE c.citationstr=?"""
            pgcur.execute(query, (self.citestr,))
            res = pgcur.fetchone()
            if res != None:
                citation_id = res[0]
        else:
            query = """SELECT citation_id
                FROM citations c
                WHERE (c.citationstr IS NULL OR c.citationstr='') AND
                    (c.url IS NULL OR c.url='') AND
                    (c.doi IS NULL OR c.doi='') AND c.authordisplay=?"""
            pgcur.execute(query, (self.authordisp,))
            res = pgcur.fetchone()
            if res != None:
                citation_id = res[0]
        
        # If a matching citation was not found, create one in the database.
        if citation_id == None:
            query = """INSERT INTO citations
                (citationstr, url, doi, authordisplay)
                VALUES (?, ?, ?, ?)"""
            pgcur.execute(query, (self.citestr, self.url, self.doi, self.authordisp))
            citation_id = pgcur.lastrowid
            if do_commit:
                pgcur.connection.commit()

        return citation_id

