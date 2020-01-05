
try:
    from ConfigParser import RawConfigParser #python2
except:
    from configparser import RawConfigParser #python3
import os.path as path
import re


class ConfigError(Exception):
    """
    A basic exception class for reporting configuration errors in the taxonomy
    configuration file.
    """
    def __init__(self, msg):
        msg = 'Configuration file error:\n  ' + msg
        Exception.__init__(self, msg)


class ColumnFilter:
    """
    Implements filtering of a CSV row based on the value of a particular column.
    """
    def __init__(self, colname, valuelist):
        """
        The argument colname specifies the name of the CSV file column to match;
        value list specifies the list of values to try to match.  If a value in
        valuelist is prefixed by "regex:", then it is interpreted as a regular
        expression.  Otherwise, the values are interpreted as plain strings that must
        match the column value exactly.
        """
        self.column = colname
        
        # This is a hack to get the type of compiled regex objects.  The Python documentation
        # says these should be of type re.RegexObject, but no such class is defined in the
        # re module.  I consider this to be a problem with the Python implementation, so this
        # stupid hack is currently the safest way to get the object type.
        self.reobjtype = type(re.compile(''))

        self.filtervals = self._processValueList(valuelist)

    def _processValueList(self, vallist):
        # Process each search string.
        for index in range(len(vallist)):
            # Remove leading and trailing whitespace.
            vallist[index] = vallist[index].strip()

            # Check each search string to see if it is a regular expression.
            if vallist[index][0:6] == 'regex:':
                restr = vallist[index][6:]
                vallist[index] = re.compile(restr)

        return vallist

    def checkMatch(self, csvrow):
        """
        Check if a row from a CSV file matches the filter.  Returns TRUE if there is a
        match, FALSE otherwise.
        """
        testval = csvrow[self.column]

        for val in self.filtervals:
            if isinstance(val, self.reobjtype):
                if val.search(testval) != None:
                    return True
            else:
                if val == testval:
                    return True

        return False


class TaxonomyConfig(RawConfigParser):
    """
    Parses taxonomy configuration files and provides convenient access to the
    configuration values.  Several methods are provided that bundle related
    configuration values and return them as a tuple, but all configuration values
    can also be accessed using the usual ConfigParser methods (e.g., get(),
    getint(), etc.).
    """
    def optionxform(self, val):
        """
        Defines an optionxform() method that doesn't change the case of option names.
        (The default implementation makes all option names lower case.)
        """
        return str(val)
    
    def read(self, filename):
        """
        Reads the configuration information and checks that the configuration file
        contains all required information for processing CSV taxonomy files.
        """
        # Call the superclass read() method.
        filesread = RawConfigParser.read(self, filename)

        if len(filesread) == 0:
            raise ConfigError('The configuration file ' + filename + ' could not be opened.')
        else:
            self.checkConfig()

        self.confdir = path.dirname(path.abspath(filename))

    def checkConfig(self):
        """
        Performs some basic checks to make sure the configuration file is valid.  If
        any problems are found, a ConfigError exception is thrown.
        """
        # Main configuration options.
        if not(self.has_option('main', 'inputcsv')):
            raise ConfigError('No input CSV file specified.')
        if not(self.has_option('main', 'ranksys')):
            raise ConfigError('No taxonomic rank system specified.')

        # Taxonomy configuration options.
        if not(self.has_option('taxonomy', 'taxonomyid')):
            raise ConfigError('No taxonomy ID specified.')
        if not(self.has_option('taxonomy', 'name')):
            raise ConfigError('The taxonomy name is missing.')
        if not(self.has_option('taxonomy', 'ismaster')):
            raise ConfigError('Master taxonomy status not specified.')
    
        # Root taxon information.
        if not(self.has_option('root', 'rank')):
            raise ConfigError('Missing root taxon rank.')
        if not(self.has_option('root', 'name')):
            raise ConfigError('Missing root taxon name string.')
    
        # Citation information.
        if not(self.has_option('citation', 'citationstr')):
            raise ConfigError('The full citation string is missing.')
        if not(self.has_option('citation', 'authordisplay')):
            raise ConfigError('The short citation string (authordisplay) is missing.')
    
        # Rank name to CSV column name mappings.
        if not(self.has_section('rank_mappings') and len(self.options('rank_mappings')) > 0):
            raise ConfigError('The "rank_mappings" section is missing or empty.')

    def getMainSettings(self):
        """
        Returns the top-level configuration information.  If no encoding was specified in
        the configuration file, 'utf-8' is returned by default.
        """
        ranksys = self.getint('main', 'ranksys')

        # Build an absolute path to the CSV file, using the location of the configuration
        # file as the starting point for relative paths.
        inputcsv = self.get('main', 'inputcsv')
        if not(path.isabs(inputcsv)):
            inputcsv = path.abspath(path.join(self.confdir, inputcsv))

        if self.has_option('main', 'encoding'):
            encoding = self.get('main', 'encoding')
        else:
            encoding = 'utf-8'
        return (ranksys, inputcsv, encoding)

    def getAcceptFilter(self):
        """
        Returns a tuple containing the name of the filter column and a list of values which
        indicate a taxon should be included in the taxonomy.  If this section is missing, or
        if either of the filter settings are missing, or if no filter column name is provided,
        then (None, None) is returned.
        """
        return self._getFilter('acceptvalues')

    def getRejectFilter(self):
        """
        Returns a tuple containing the name of the filter column and a list of values which
        indicate a taxon should be included in the taxonomy.  If this section is missing, or
        if either of the filter settings are missing, or if no filter column name is provided,
        then (None, None) is returned.
        """
        return self._getFilter('rejectvalues')

    def _getFilter(self, optname):
        """
        Internal method for building a list of accept or reject column search values.  This method
        expects that the required settings are in a section called "taxafilter" with the filter
        column specified by the option "filtercolumn" and the filter values option specified by the
        value of the optname argument.
        """
        filterobj = None

        if self.has_section('taxafilter'):
            if self.has_option('taxafilter', 'filtercolumn') and self.has_option('taxafilter', optname):
                filtercol = self.get('taxafilter', 'filtercolumn')

                accstrs = self.get('taxafilter', optname)
                accstrs = accstrs.split(',')

                if filtercol != '':
                    filterobj = ColumnFilter(filtercol, accstrs)

        return filterobj

    def getTransformations(self):
        """
        Returns a dictionary mapping CSV column names to regular expressions that should
        be used to transform the values of the column strings.
        """
        transforms = {}

        if self.has_section('transformations'):
            for namepair in self.items('transformations'):
                restr = namepair[1].strip('"')
                transforms[namepair[0]] = re.compile(restr)

        return transforms

    def getTaxonomySettings(self):
        """
        Returns all taxonomy-specific settings.
        """
        taxonomyid = self.getint('taxonomy', 'taxonomyid')
        taxonomyname = self.get('taxonomy', 'name')
        ismaster = self.getboolean('taxonomy', 'ismaster')

        return (taxonomyid, taxonomyname, ismaster)

    def getRootSettings(self):
        """
        Returns the settings that specify the root taxon_concept for the taxonomy.
        """
        rootrank = self.get('root', 'rank')
        rootname = self.get('root', 'name')
        roottaxoid = self.get('root', 'taxonomyid')

        return (rootrank, rootname, roottaxoid)

    def getCitationSettings(self):
        """
        Returns all four taxonomy citation-related settings as a tuple.  The options "url"
        and "doi" are not required, and if they are not specified, their values are returned
        as None.
        """
        citestr = self.get('citation', 'citationstr')
        citeauthor = self.get('citation', 'authordisplay')
        citedoi = None
        if self.has_option('citation', 'doi'):
            citedoi = self.get('citation', 'doi')
        citeurl = None
        if self.has_option('citation', 'url'):
            citeurl = self.get('citation', 'url')
        
        return (citestr, citeauthor, citedoi, citeurl)

    def getTaxaCitationSettings(self):
        """
        Returns citation settings for taxon names.  Returns a tuple that contains
        three items, in the following order: 1) a tuple containing, in order, the
        CSV column name that contains rank-inspecific short author strings, and the
        CSV column name that contains the full citations for the author strings; 2)
        a dictionary mapping MOL rank names to CSV column names that contain
        rank-specific citation information; and 3) a Boolean value that indicates
        whether the parser should to try to fix casing problems.  If the first
        setting is not defined, '' is returned, and if no rank-specific citation
        information is defined, an empty dictionary is returned.
        """
        nameauthcol = ''
        if self.has_option('citation', 'nameauthorcol'):
            nameauthcol = self.get('citation', 'nameauthorcol')
        namefullcitecol = ''
        if self.has_option('citation', 'namefullcitecol'):
            namefullcitecol = self.get('citation', 'namefullcitecol')

        rankciteinfo = {}
        for optname, optval in self.items('citation'):
            if optname.find('citation_') == 0:
                rankname = optname.replace('citation_', '')
                rankciteinfo[rankname] = optval

        fixcasing = False
        if self.has_option('citation', 'fixcasing'):
            fixcasing = self.getboolean('citation', 'fixcasing')

        return ((nameauthcol, namefullcitecol), rankciteinfo, fixcasing)

    def getSynonymSettings(self):
        """
        Returns the two synonym settings as a tuple.  If a separator character is not defined,
        ',' is used.  If the synonym settings are not defined, (None, None) is returned.
        """
        colname = sep = None

        if self.has_section('synonyms'):
            if self.has_option('synonyms', 'colname'):
                colname = self.get('synonyms', 'colname')
                if self.has_option('synonyms', 'separator'):
                    sep = self.get('synonyms', 'separator')
                else:
                    sep = ','

        return (colname, sep)

    def getRankMappings(self, rankt, ranksys):
        """
        Returns a dictionary mapping rank IDs in the taxonomy database to column names
        in the CSV file.  The rank IDs are used as the dictionary keys, and the CSV file
        column names are the values.  An appropriate RankTable must be provided.
        """
        ranktoCSV = {}
        for namepair in self.items('rank_mappings'):
            try:
                rankID = rankt.getID(namepair[1], ranksys)
            except KeyError as e:
                raise ConfigError('"' + str(e.args[0]) + '" is not a valid taxonomic rank name.')

            ranktoCSV[rankID] = namepair[0]
        #print ranktoCSV

        return ranktoCSV
        
