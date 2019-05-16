
import csv
from taxacomponents import RankTable, Taxon, Name
from nameresolve import NamesResolver


class TaxoCSVError(Exception):
    """
    A basic exception class for reporting errors encountered while processing a
    taxonomy CSV file.
    """
    def __init__(self, msg):
        msg = 'Error while parsing taxonomy CSV file:\n  ' + msg
        Exception.__init__(self, msg)


class UnicodeDictReader(csv.DictReader):
    """
    Extends the DictReader class to support various encodings of input CSV files.  All
    data are parsed into Python unicode strings.  By default, data are assumed to be
    encoded using UTF-8, but alternative encodings can also be specified.
    """
    def __init__(self, csvfile, fieldnames=None, restkey=None, restval=None, dialect='excel',
            encoding='utf-8', *args, **kwds):
        """
        This initializer has the same signature as that of standard csv.DictReader object,
        except that an addiitonal argument is added to specify the source encoding.
        """
        self.encoding = encoding

        # Call the superclass constructor.
        csv.DictReader.__init__(self, csvfile, fieldnames, restkey, restval, dialect,
                *args, **kwds)

    def next(self):
        row = csv.DictReader.next(self)
        for key in row.keys():
            row[key] = unicode(row[key], self.encoding)

        return row


class UnicodeDictWriter(csv.DictWriter):
    """
    Extends the DictWriter class to support unicode strings.  All output strings are
    encoded as UTF-8 text.
    """
    def writerow(self, row):
        encrow = {}
        for key in row.keys():
            if isinstance(row[key], basestring):
                # Only try to encode instances of string types.
                encrow[key] = row[key].encode('utf-8')
            else:
                encrow[key] = row[key]

        # Call the superclass method.
        csv.DictWriter.writerow(self, encrow)


class CSVTaxonomyParser:
    """
    Reads taxonomy information from a CSV file where each row represents a single
    taxonomic unit along with its higher taxonomy.  Builds a tree of Taxon objects
    that represents the taxonomy contained in the CSV file.
    """
    def __init__(self):
        self.totalrows = 0
        self.totaltaxa = 0

    def parseCSV(self, taxoconfig, dbcur):
        """
        Parses a taxonomy from a CSV file.  Requires a valid TaxonomyConfig object and
        database cursor.  The configuration object is required to provide the location
        of the CSV file to parse along with the metadata needed to properly build the
        taxon tree.  The database cursor is needed to access taxonomic rank information.
        Returns the root taxon for the taxon tree structure.
        """
        self.tc = taxoconfig

        # Get the top-level configuration settings.
        ranksys, inputcsv, charencoding = self.tc.getMainSettings()
        
        # Get the taxonomy ID from the configuration file.
        taxonomyid = self.tc.getTaxonomySettings()[0]

        # Get the transformations dictionary.
        self.transforms = self.tc.getTransformations()
        
        # Initialize the rank lookup table.
        rankt = RankTable()
        rankt.loadFromDB(dbcur)

        # Get the synonyms settings.
        self.syn_col, self.syn_sep = self.tc.getSynonymSettings()

        # Get the taxon name citation information.
        self.nameciteinfo = self.tc.getTaxaCitationSettings()

        # Instantiate a NamesResolver to use for processing name strings and author
        # display strings.
        self.resolver = NamesResolver()

        # Get the CSV file column name mappings from the configuration file.
        ranktoCSV = self.tc.getRankMappings(rankt, ranksys)
        
        # Create a Taxon object for the root of the taxonomy.
        rootrank, rootname, roottaxo_id = self.tc.getRootSettings()
        taxonomyroot = Taxon(taxonomyid, rankt.getID(rootrank, ranksys), rankt, 0, rootname, roottaxo_id, True)
        #print taxonomyroot

        # Open and parse the input CSV file.
        with open(inputcsv, 'rU') as fin:
            reader = UnicodeDictReader(fin, encoding=charencoding)
            self._readCSVRows(reader, taxonomyroot, ranktoCSV, rankt)

        return taxonomyroot

    def _extractGenus(self, spstr):
        """
        Parses out and returns the genus name from a scientific name string.
        """
        pieces = spstr.strip().split(' ')

        return pieces[0]

    def _extractSpecies(self, spstr):
        """
        Parses out and returns the species name from a scientific name string.
        If the string does not contain a species name, '' is returned.
        """
        pieces = spstr.strip().split(' ')

        if len(pieces) < 2:
            return ''
        else:
            return pieces[0] + ' ' + pieces[1]

    def _extractSubsp(self, spstr):
        """
        Parses out and returns the subspecies name from a scientific name string.
        If the string does not contain a subspecies name, '' is returned.
        """
        pieces = spstr.strip().split(' ')

        if len(pieces) < 3:
            return ''
        else:
            return pieces[0] + ' ' + pieces[1] + ' ' + pieces[2]

    def _applyTransforms(self, row):
        """
        Applies any column value transformations that were defined in the configuration file.
        """
        for colname, regex in self.transforms.iteritems():
            row[colname] = regex.sub('', row[colname])

    def _addSynonyms(self, taxon, row):
        """
        Parses synonymous names from the CSV row.  Returns a list of Name objects,
        one for each synonym.
        """
        if self.syn_col != None:
            rawvals = row[self.syn_col].split(self.syn_sep)

            # Verify that we actually got synonym data.
            if len(rawvals) == 1 and rawvals[0] == '':
                rawvals = []

            for val in rawvals:
                name = Name(val.strip())
                taxon.namelist.addName(name, False, False)

    def _addNameCitation(self, taxon, authorstr, fullcitestr):
        """
        Attaches citation information to the name of a taxon.  The argument authorstr
        provides the author display string for the name authority citation; fullcitestr
        provides the full citation string.
        """
        fullcitestr = fullcitestr.strip()

        if authorstr.strip() != '':
            authorinfo = self.resolver.processAuthorString(authorstr)

            if self.nameciteinfo[2]:
                # Attempt to fix the casing of the author string.
                authorinfo['authorstr'] = authorinfo['authorstr'].title()

            taxon.name.updateCitation(fullcitestr, authorinfo['authorstr'])
            taxon.setUseParens(authorinfo['hasparens'])
        elif fullcitestr != '':
            taxon.name.updateCitation(fullcitestr, '')

    def _readCSVRows(self, reader, taxonomyroot, ranktoCSV, rankt):
        """
        Reads taxonomy information from a CSV file and attaches the taxonomic unit
        for each row as a descendent of taxonomyroot.
        """
        self.totalrows = self.totaltaxa = 0

        # Get the column accept/reject filters.
        acceptfilter = self.tc.getAcceptFilter()
        rejectfilter = self.tc.getRejectFilter()

        # Get the name of the column containing the species name string so we can use it to
        # parse out the genus name, if needed.
        ranksys = self.tc.getMainSettings()[0]
        sp_rankid = rankt.getID('Species', ranksys)
        sp_colname = ranktoCSV[sp_rankid]

        # See if we also are parsing subspecies names from the species name string.  If so,
        # we need to preprocess species name strings to remove the subspecies name.
        parse_subsp = False
        ssp_rankid = rankt.getID('Subspecies', ranksys)
        if ssp_rankid in ranktoCSV:
            parse_subsp = ranktoCSV[ssp_rankid] == '_ParseSubspecies'

        # Process each row in the CSV file, building up the tree of Taxon objects as we go.
        rankorder = sorted(ranktoCSV.keys())
        for row in reader:
            self.totalrows += 1

            # If a filter column and values were specified, check them to see if this row should
            # be included in the taxonomy.
            if acceptfilter != None:
                if not(acceptfilter.checkMatch(row)):
                    # Skip to the next row.
                    continue
            if rejectfilter != None:
                if rejectfilter.checkMatch(row):
                    # Skip to the next row.
                    continue

            # Apply any column value transformations.
            self._applyTransforms(row)

            # Start with the root as the parent.
            curparent = taxonomyroot
        
            # Process each rank in hierarchical order.
            for rankid in rankorder:
                # Get the name string.
                try:
                    if ranktoCSV[rankid] == '_ParseGenus':
                        # Parse out the genus name from the species name string.
                        namestr = self._extractGenus(row[sp_colname])
                    elif ranktoCSV[rankid] == '_ParseSubspecies':
                        # Parse out the subspecies name from the species name string.
                        namestr = self._extractSubsp(row[sp_colname])
                    else:
                        namestr = row[ranktoCSV[rankid]].strip()
                        if (rankid == sp_rankid) and parse_subsp:
                            namestr = self._extractSpecies(namestr)
                except KeyError as e:
                    raise TaxoCSVError('The column "' + e.args[0] + '" was not found in the taxonomy CSV file.')

                if namestr != '':
                    # Clean up the name string.
                    namestr = self.resolver.cleanNameString(namestr)
        
                    # See if this taxon has already been processed.  If not, create a Taxon object for it.
                    childtaxon = curparent.findChild(rankid, namestr)
                    if childtaxon == None:
                        childtaxon = curparent.createChild(rankid, namestr)
                        # If we have rank-specific citation data for this taxon, attach it.
                        rankname = rankt.getName(rankid)
                        if rankname in self.nameciteinfo[1]:
                            self._addNameCitation(childtaxon, row[self.nameciteinfo[1][rankname]], '')
                        self.totaltaxa += 1
        
                    # Set the child Taxon as the parent for the next taxon to process in the row.
                    curparent = childtaxon

            # Attach any synonyms to the last child taxon we processed for this row.
            if self.syn_col != None:
                syns = self._addSynonyms(childtaxon, row)

            # Attach any rank-inspecific citation data to the name of the last child taxon we
            # processed for this row.
            nameauthstr = ''
            namefullcitestr = ''
            if self.nameciteinfo[0][0] != '':
                nameauthstr = row[self.nameciteinfo[0][0]]
            if self.nameciteinfo[0][1] != '':
                namefullcitestr = row[self.nameciteinfo[0][1]]
            if (nameauthstr != '') or (namefullcitestr != ''):
                self._addNameCitation(childtaxon, nameauthstr, namefullcitestr)


    def getStats(self):
        """
        Returns the total number of CSV file rows and unique taxa that were processed.
        The values are returned as a tuple.
        """
        return (self.totalrows, self.totaltaxa)

