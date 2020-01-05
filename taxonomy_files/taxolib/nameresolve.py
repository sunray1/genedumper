
try:
    import urllib, urllib.request, urllib.error
except:
    import urllib, urllib2
try:
    from BaseHTTPServer import BaseHTTPRequestHandler
except:
    from http.server import BaseHTTPRequestHandler
import socket
import json
import xml.etree.ElementTree as et
import re
import sys, time
import pprint
from taxolib.taxacomponents import Taxon, Citation
from taxolib.taxonvisitor import TaxonVisitor


def getResolversList():
    """
    A simple factory function that instantiates each concrete resolver class
    and returns the resolver objects in a list.
    """
    resolvers = []
    resolvers.append(CoLNamesResolver())
    resolvers.append(ZoobankNamesResolver())

    return resolvers


class NamesResolver(TaxonVisitor):
    def __init__(self, numtaxa=-1, maxdepth=-1):
        # Call the superclass initializer.
        TaxonVisitor.__init__(self, numtaxa, maxdepth)

        # The initial timeout limit (in seconds) for a request.
        self.timeout = 60

        # If an initial request fails due to a connection timeout, subsequent
        # retries can use progressively increasing timeout limits.  Each retry,
        # the timeout limit will by multiplied by self.timeoutfactor.
        self.timeoutfactor = 1.5

        # The number of times to retry a request in case of TCP connection failure
        # or an HTTP error.
        self.maxretries = 2

        # The delay (in seconds) between retry attempts.
        self.retrydelay = 20

        # By default, do not retry a request after getting an HTTP 404 error.
        self.retry404 = False

        # Get the HTTP response code messages.
        self.HTTPresponses = BaseHTTPRequestHandler.responses

        # A regular expression to match runs of 2 or more whitespace characters.
        self.wsregx = re.compile('\s{2,}')

        # Options to control whether only names that do not already exist in the
        # database with authorship citation data are resolved.  For this functionality,
        # a comparison taxonomy must be provided for name matching.
        self.comparison_taxonomy = None
        self.skip_if_existing = False

    def getSourceDescription(self):
        """
        Return a short text string describing the data source for this NameResolver.
        This method should be overridden by child classes.
        """
        return ''

    def setComparisonTaxonomy(self, taxonomy_id):
        """
        Sets a comparison taxonomy to use for checking if a name already exists in
        the taxonomy database.
        """
        self.comparison_taxonomy = taxonomy_id

    def setSkipIfExisting(self, skip):
        """
        If skip is True, then the name resolver will check if the target name already
        exists in the database as part of a comparison taxonomy.  If so, and if it
        already has name citation data, then no further attempt will be made to look
        up citation data for the name.
        """
        self.skip_if_existing = skip

    def _checkExistingCiteData(self, taxon):
        """
        This method checks whether a name already exists in the database as part of a
        comparison taxonomy, and if so, whether it has citation data.  A search is
        performed only if self.skip_if_existing is True.  If a match is found, this
        method will load the name/citation data into the target taxon concept and
        specify whether parentheses should be used when displaying the author display
        string.  If this is all successful, then the method returns True; otherwise,
        it returns False.
        """
        if self.skip_if_existing and self.comparison_taxonomy != None:
            res = Taxon.find(self.pgcur, taxon.name.namestr, taxon.rankt,
                    taxonomy_id=self.comparison_taxonomy, rank_id=taxon.rank_id, pref_names_only=True)
            if len(res) == 1:
                ctaxon = res[0]
                if ctaxon.name.getCitation() != None:
                    taxon.name.loadFromDB(self.pgcur, ctaxon.name.idnum)
                    taxon.setUseParens(ctaxon.getUseParens())
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

    def queryJSON(self, queryurl):
        """
        Queries a Web service that returns results in JSON format.  If the
        query is successful, the results are parsed as JSON and returned as a
        Python object.  See the documentation for the method _HTTPQuery() for
        more information about how TCP and HTTP errors are handled.
    
        Arguments:
          url -- The base URL of the Web service.
          querystr -- The query string to use.
    
        Returns:  A Python object that represents the JSON query result.
        """
        res = self._HTTPQuery(queryurl)
        jsonstr = res.read()

        # Zoobank sometimes returns invalid JSON.  So far, I have seen JSON string tokens
        # with literal tab characters, which is invalid according to the standard.  Try
        # to fix some of these problems.
        jsonstr = jsonstr.replace('\t', '\\t')

        return json.loads(jsonstr)

    def queryXML(self, queryurl):
        """
        Queries a Web service that returns results in XML format.  If the
        query is successful, the results are parsed as XML and returned as a
        Python object.  XML is parsed using the Python ElementTree API.  See
        the documentation for the method _HTTPQuery() for more information
        about how TCP and HTTP errors are handled.
    
        Arguments:
          url -- The base URL of the Web service.
          querystr -- The query string to use.
    
        Returns:  An ElementTree object that represents the XML data tree.
        """
        res = self._HTTPQuery(queryurl)

        return et.parse(res)

    def _HTTPQuery(self, queryurl):
        """
        Queries a Web service.  If the query is successful, a reference to the
        results is returned.  If TCP or HTTP errors are encountered, this method
        will retry the request after waiting self.retrydelay seconds.  If the
        total number of retries exceeds self.maxretries and the request is still
        unsuccessful, the TCP or HTTP exception will be passed up the call chain.
        Because Python is inconsistent in how it reports timeout errors (they are
        sometimes raised as a bare socket.timeout object, and sometimes come as a
        socket.timeout wrapped in a URLError object), this method detects both
        ways of reporting timeouts and will raise them only as socket.timeout
        objects so that client code doesn't have to worry about inconsistent
        timeout error reporting.
    
        Arguments:
          url -- The base URL of the Web service.
          querystr -- The query string to use.
    
        Returns:  A reference to the request result.
        """
        # Initialize variables for managing request retry attempts.
        success = False
        retries = 0
        timeoutfailures = 0

        # Try the request until we get a result back or the maximum number of
        # retries has been exceeded due to TCP or HTTP errors.
        while not(success):
            try:
                res = urllib2.urlopen(queryurl, None, self.timeout * (self.timeoutfactor**timeoutfailures))
                success = True
            except (socket.error, urllib2.URLError) as err:
                # The exception handling code is complex because there are a lot of things
                # that can go wrong, and worse, because Python is inconsistent in how it
                # reports these errors.  For example, timeout errors sometimes come as a
                # bare socket.error object, and sometimes they come as a socket.error wrapped
                # in a urllib2.URLError object.

                # First, check if we got a timeout error.
                gottimeout = ((type(err) is socket.timeout) or
                        ((type(err) is urllib2.URLError) and (type(err.reason) is socket.timeout)))

                retries += 1
                if retries > self.maxretries:
                    print('Connection attempt failed; maximum retries exceeded.')
                    # Pass the exception on to the caller.  Handle timeout errors by only
                    # passing on a socket.timeout object so that client code doesn't have
                    # to deal with multiple timeout error types.
                    if gottimeout and (type(err) is urllib2.URLError):
                        raise err.reason
                    else:
                        raise

                # Determine what sort of error occured, print an appropriate error message
                # and, if needed, adjust the timeout limit.
                if type(err) is urllib2.HTTPError:
                    if err.code == 404 and not(self.retry404):
                        raise
                    else:
                        print('HTTP error ' + str(err.code) + ': ' + self.HTTPresponses[err.code][0] + '.')
                elif gottimeout:
                    print('Connection error: connection attempt timed out.')
                    # Keep track of the number of timeout failures so that the timeout limit
                    # can be properly adjusted.
                    timeoutfailures += 1
                    print ('Increasing socket timeout limit to ' +
                            str(self.timeout * (self.timeoutfactor**timeoutfailures)) + ' seconds.')
                elif type(err) is socket.error:
                    print('TCP socket error:', err)
                else:
                    print('Connection error:', err.reason)

                print('Retrying request in ' + str(self.retrydelay) + ' seconds.')
                time.sleep(self.retrydelay)
    
        return res

    def resolve(self, pgcur, taxon):
        """
        This method merely saves a database cursor object for use by other methods of the
        class during name resolution and then calls the visit() method of TaxonVisitor.
        """
        self.pgcur = pgcur

        self.visit(taxon)

    def processAuthorString(self, authorstr):
        """
        Standardizes author display strings so that they are more easily comparable
        across citation data sources.  Also checks for wrapping parentheses and removes
        them, if found.  Returns a dictionary with two elements: 1) 'authorstr', which
        contains the author display string; and 2) 'hasparens', a boolean value that
        indicates whether wrapping parentheses were used.
        """
        authorstr = self._commonStringCleanup(authorstr)

        # Remove commas.
        authorstr = authorstr.replace(',', '')

        # Replace the "ae" grapheme (e.g., sometimes used in "Linnaeus") with the
        # letters "ae".
        authorstr = authorstr.replace(u'\xe6', 'ae')

        # Check for and remove wrapping parentheses.
        hasparens = (authorstr != '' and authorstr[0] == '(' and authorstr[-1] == ')')
        if hasparens:
            # Remove the parentheses from the author string.
            authorstr = authorstr[1:-1]

        return { 'authorstr': authorstr, 'hasparens': hasparens }
    
    def cleanNameString(self, namestr):
        """
        Standardizes taxon name strings so that they are more easily comparable across
        citation data sources.
        """
        namestr = self._commonStringCleanup(namestr)

        # Make sure only the first letter of the name string is capitalized.
        namestr = namestr[0].upper() + namestr[1:].lower()

        # If this name includes a subgenus designation, make sure the first character of
        # the subgenus name is capitalized.
        s_index = namestr.find('(')
        e_index = namestr.find(')')
        if (s_index != -1) and (e_index != -1):
            namestr = namestr[:s_index+1] + namestr[s_index+1].upper() + namestr[s_index+2:]

        return namestr

    def cleanCiteString(self, citestr):
        """
        Standardizes citation strings so that they are more easily comparable across
        citation data sources.
        """
        citestr = self._commonStringCleanup(citestr)

        # Replace the en dash with a hyphen.
        fullcite = citestr.replace(u'\u2013', '-')

        return citestr

    def _getSpeciesSearchStrs(self, namestr):
        """
        Returns a list of name strings to try when searching for a matching taxon in
        a resolver resource.  If namestr does not include a subgenus, only namestr
        is returned.  Otherwise, the list will contain, in order: the original name
        (namestr), the name without the subgenus, the name without the genus (that is,
        considering the subgenus to be the genus).  This method assumes that the
        components of namestr are separated by only a single whitespace character;
        that is, that they have been processed by _commonStringCleanup().
        """
        nameslist = [namestr]

        s_index = namestr.find('(')
        e_index = namestr.find(')')
        if (s_index != -1) and (e_index != -1):
            subgenus = namestr[s_index+1 : e_index]
            genus = namestr[:s_index-1]
            sppart = namestr[e_index+2:]
            nameslist.append(genus + ' ' + sppart)
            nameslist.append(subgenus + ' ' + sppart)

        return nameslist

    def _commonStringCleanup(self, datastr):
        """
        Common string standardization code used by the clean*String() methods.
        """
        # Remove leading/trailing whitespace.
        datastr = datastr.strip()

        # Replace runs of 2 or more whitespace characters with a single space.
        datastr = self.wsregx.sub(' ', datastr)

        # Remove spaces before commas.
        datastr = datastr.replace(' ,', ',')

        return datastr


class ZoobankNamesResolver(NamesResolver):
    """
    A names resolver that attempts to retrieve name citation information from Zoobank.
    A name usage record from Zoobank is considered to match the target name if both
    the name string and taxonomic rank name match.  If short author information
    already exists for the target name, it is used to further confirm a match.  The
    initial Zoobank search is by name string only (this is all that the Zoobank API
    supports).  Unfortunately, checking a large number of names against Zoobank is
    very slow because the Zoobank API only supports searching for a single name or
    UUID at one time and each query is slow.
    """
    def __init__(self, numtaxa=-1, maxdepth=-1):
        # Call the superclass initializer.
        NamesResolver.__init__(self, numtaxa, maxdepth)

        self.zoobank_name_url = 'http://zoobank.org/NomenclaturalActs.json/'
        self.zoobank_ref_url = 'http://zoobank.org/References.json/'

    def getSourceDescription(self):
        return 'Zoobank'

    def resolve(self, pgcur, taxon):
        # Initialize state-tracking variables.
        self.last_lookup_rank = ''
        self.last_lookup_succeeded = False

        # Some reference UUIDs provided by Zoobank cause queries to completely hang
        # and eventually fail.  Keep track of failed UUIDs to avoid wasting time
        # repeatedly attempting to retrieve them.
        self.failedrefIDs = []

        # Call the superclass implementation.
        NamesResolver.resolve(self, pgcur, taxon)

    def doRecursion(self):
        """
        If a genus name lookup to Zoobank fails to return any results, there is no
        reason to waste time looking up the species names within the genus.  This
        method checks for this condition and prevents the tree recursion from
        continuing to the child taxa if the genus name lookup failed.
        """
        if self.last_lookup_rank == 'Genus' and not(self.last_lookup_succeeded):
            return False
        else:
            return True

    def searchZoobankForTaxon(self, taxon, name_searchstr):
        """
        Searches Zoobank for a taxon name string.  If a match is found, returns a
        tuple containing the parsed JSON result, the cleaned name string, author
        information, and the full citation string; otherwise, returns None.
        """
        name = taxon.name

        # The Zoobank API returns 404 if a search returns no results, so we need to
        # check for this.  I've also found that the service is somewhat unreliable.
        # Requests will occasionally time out, even after multiple retries and with
        # rather long timeout limits (e.g., 90 seconds).  Rather than let these
        # occasional failures crash the program, catch these timeout errors, print a
        # message to the console, and move on.
        queryurl = self.zoobank_name_url + urllib.quote(name_searchstr.replace(' ', '_'))
        try:
            rjson = self.queryJSON(queryurl)
        except urllib2.HTTPError as e:
            if e.code == 404:
                rjson = []
        except ValueError as e:
            print ('Attempt to search Zoobank for the name "' + name_searchstr
                    + '" failed because the result was not a valid JSON string.')
            rjson = []
        except socket.timeout as e:
            print ('Attempt to search Zoobank for the name "' + name_searchstr
                    + '" failed due to multiple timeout errors.')
            rjson = []
        #print rjson
        
        match = None
        matchcnt = 0

        # Inspect each result and see if we have a match.
        for res in rjson:
            # Split the Zoobank "protonym" into the name string and author string.
            #print res['cleanprotonym']
            #print res['namestring']
            #print res['rankgroup']
            parts = res['cleanprotonym'].rsplit(res['namestring'], 1)
            #print parts

            # Verify that the name string was found in the protonym string.  If not,
            # then either one or both was misspelled and we can't use them.
            if len(parts) == 2:
                authorinfo = self.processAuthorString(parts[1])
                rnamestr = self.cleanNameString(parts[0] + res['namestring'])
                #print authorinfo
                #print rnamestr

                # See if we have a match.
                havematch = False
                if res['rankgroup'] == taxon.getRankString() and rnamestr == name_searchstr:
                    if name.getCitation() != None:
                        # If the existing citation already has short author information, use it
                        # to further verify the match.
                        if name.getCitation().authordisp == authorinfo['authorstr']:
                            havematch = True
                    else:
                        havematch = True

                if havematch:
                    # Keep track of the total number of name/rank matches.  If there is more than
                    # one, we don't know which is the correct match and so can't use any.
                    matchcnt += 1
                    match = (res, rnamestr, authorinfo)

        if matchcnt == 1:
            authorinfo = match[2]
            authorprint = authorinfo['authorstr']
            if authorinfo['hasparens']:
                authorprint = '(' + authorprint + ')'
            print('  Found match:', match[1], authorprint)

            # See if Zoobank contains detailed publication information.
            # Frustratingly, Zoobank stores publication DOIs but does not return them as
            # part of their API results.  We can get the full citation string, though.
            fullcite = ''
            if (match[0]['OriginalReferenceUUID'] != '' and
                    match[0]['OriginalReferenceUUID'] not in self.failedrefIDs):
                # I have found that some reference requests to the Zoobank API (using supposedly
                # legitimate UUIDs from Zoobank) never return.  Thus, we need to watch for
                # TCP timeout failures here and not let them terminate the program.
                try:
                    queryurl = self.zoobank_ref_url + urllib.quote(match[0]['OriginalReferenceUUID'])
                    rjson = self.queryJSON(queryurl)
                    #pprint.pprint(rjson)
                    # Get the full citation string.
                    fullcite = self.cleanCiteString(rjson[0]['value'])
                    print(fullcite)
                except socket.timeout as e:
                    # Keep track of failed IDs so we don't waste time on them later.
                    self.failedrefIDs.append(match[0]['OriginalReferenceUUID'])
                    print ('Attempt to retrieve Zoobank reference UUID ' +
                            match[0]['OriginalReferenceUUID'] + ' failed due to multiple timeout errors.')

            return match + (fullcite,)
        else:
            return None

    def processTaxon(self, taxon, depth):
        """
        Uses Zoobank to try to retrieve author and citation information for taxa
        name strings.
        """
        # Update the rank state-tracking variable.
        self.last_lookup_rank = taxon.getRankString()

        # Only operate on taxa of family-group rank or lower.  For higher taxa,
        # Zoobank only stores "Higher" as the rank name, so we cannot make rank
        # comparisons.
        if taxon.rank_id < taxon.rankt.getID('Family', taxon.ranksys):
            return

        # See if we should skip this taxon because its name already exists in the
        # database and has citation data.
        if self._checkExistingCiteData(taxon):
            # Consider this a lookup "success" so that we recurse to the next level,
            # even though lookup was skipped for this taxon.
            self.last_lookup_succeeded = True
            print('Skipping "' + taxon.name.namestr + '"; citation data already in database.')
            return

        print('Attempting to resolve:', taxon.name.namestr)
        sys.stdout.flush()

        # Search for this species' name in Zoobank, including name variants if it has a
        # subgenus designation.
        if taxon.rank_id >= taxon.rankt.getID('Species', taxon.ranksys):
            name_searchstrs = self._getSpeciesSearchStrs(taxon.name.namestr)
            for name_searchstr in name_searchstrs:
                searchres = self.searchZoobankForTaxon(taxon, name_searchstr)
                if searchres != None:
                    break
        else:
            searchres = self.searchZoobankForTaxon(taxon, taxon.name.namestr)

        # Update the success state-tracking variable.
        if searchres != None:
            self.last_lookup_succeeded = True
        else:
            self.last_lookup_succeeded = False
            return

        res, rname, authorinfo, fullcite = searchres

        # Update the name's citation information.
        taxon.name.updateCitation(fullcite, authorinfo['authorstr'])
        taxon.setUseParens(authorinfo['hasparens'])

class CoLNamesResolver(NamesResolver):
    """
    A names resolver that attempts to retrieve name citation information from the
    Catalog of Life (CoL).  CoL does not report author information for any taxonomic
    ranks besides species-group names.  Thus, this resolver only operates on
    species-group taxa.  A record from CoL is considered to match the target name if
    the full name string and taxonomic rank name are identical and if the kingdom
    name matches the target kingdom.
    """
    def __init__(self, numtaxa=-1, maxdepth=-1):
        # Call the superclass initializer.
        NamesResolver.__init__(self, numtaxa, maxdepth)

        self.col_url = 'http://www.catalogueoflife.org/col/webservice?'

        # To increase the reliability of name matching, limit results to those within
        # a single kingdom.
        self.search_kingdom = 'Animalia'

    def getSourceDescription(self):
        return 'Catalog of Life'

    def searchCoLForTaxon(self, taxon, name_searchstr, return_no_author=False):
        """
        Searches for a taxon name string (the only type of search supported by CoL).
        If a match is found, returns a tuple containing the parsed XML result, the
        rank string, the cleaned name string, and the author information; otherwise,
        returns None.  If return_no_author is True, the method will return the
        search results even if no author information was found.
        """
        args = { 'name': name_searchstr, 'response': 'full', 'format': 'xml' }
        queryurl = self.col_url + urllib.urlencode(args)
        
        try:
            res = self.queryXML(queryurl)
        except socket.timeout as e:
            # Don't let occasional hard timeout errors crash long-running lookup efforts.
            # Print a message to the user and move on.
            print ('Attempt to search Catalog of Life for the name "' + taxon.name.namestr
                    + '" failed due to multiple timeout errors.')
            return None
    
        # Make sure we found a match.  Note that the CoL Web API documentation claims that searches
        # only return exact matches, but this is not true.  Name strings that contain the search
        # string are also returned, so we need to check each to see if we have an exact match.
        match_count = 0
    
        
        for result_tag in res.getroot():
            # Retrieve and process the taxon name.
            sname = self.cleanNameString(result_tag.find('name').text)
            # Retrieve the rank.
            srank = result_tag.find('rank').text
            # Get the kingdom name.
            kingdomnode = result_tag.find('./classification/taxon[1]/name')
            
            if kingdomnode == None:
                kingdomnode = result_tag.find('./accepted_name/classification/taxon[1]/name')
            if kingdomnode != None:
                skingdom = kingdomnode.text
            else:
                skingdom = ''

            # Check for a match.
            if (name_searchstr == sname) and (srank == taxon.getRankString() and
                    skingdom == self.search_kingdom):
                match = (result_tag, sname, srank)
                match_count += 1
            #print srank, sname, skingdom, match_count
        # If there is more than one name/rank match, we don't know which is the correct
        # match and so we can't use any.
        if match_count != 1:
            return None
        # Retrieve and process the author element.
        authtag = match[0].find('./author')
        if authtag != None:
            sauthor = authtag.text
            authorinfo = self.processAuthorString(sauthor)

            return match + (authorinfo,)
        elif return_no_author:
            return match + (None,)
        else:
            return None

    def processTaxon(self, taxon, depth):
        """
        Uses Catalog of Life to try to retrieve author information for taxa name strings.
        """
        # Only operate on species-group taxa.
        if taxon.rank_id < taxon.rankt.getID('Species', taxon.ranksys):
            return

        # See if we should skip this taxon because its name already exists in the
        # database and has citation data.
        if self._checkExistingCiteData(taxon):
            print('Skipping "' + taxon.name.namestr + '"; citation data already in database.')
            return

        print('Attempting to resolve:', taxon.name.namestr)

        # Search for this species' name in CoL, including name variants if it has a
        # subgenus designation.
        name_searchstrs = self._getSpeciesSearchStrs(taxon.name.namestr)
        for name_searchstr in name_searchstrs:
            searchres = self.searchCoLForTaxon(taxon, name_searchstr)
            if searchres != None:
                break

        if searchres == None:
            return

        res, sname, srank, authorinfo = searchres

        #print srank, sname, authorinfo

        # Update the author information for the taxon's name.
        authorprint = authorinfo['authorstr']
        if authorinfo['hasparens']:
            authorprint = '(' + authorprint + ')'
        print('  Found match:', sname, authorprint)
        taxon.name.updateCitation(None, authorinfo['authorstr'])
        taxon.setUseParens(authorinfo['hasparens'])

