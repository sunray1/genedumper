#!/usr/bin/python

import sys
from taxolib import taxodatabase
from taxolib.taxacomponents import Citation
from taxolib.taxonomy import Taxonomy
from taxolib.taxoconfig import TaxonomyConfig, ConfigError
from taxolib.csvtaxonomy import CSVTaxonomyParser, TaxoCSVError
import taxolib.nameresolve as nameresolve
from argparse import ArgumentParser


# Generate the help message for the resolvers option.
resolvers = nameresolve.getResolversList()
resolvernames = [resolver.getSourceDescription() for resolver in resolvers]
for cnt in range(len(resolvernames)):
    resolvernames[cnt] = '"' + str(cnt) + '" = ' + resolvernames[cnt]
citehelpstr = ('the name citation resolver(s) to use (' + ', '.join(resolvernames) +
    ', "all" = all resolvers [default], "none" = no citation resolution)')

argp = ArgumentParser(description='Loads taxonomies from CSV files into the MOL taxonomy database \
schema.  The single required argument provides the location of a configuration file that specifies \
the input CSV file and how it should be parsed.  The program also needs to know how to connect to \
the taxonomy database, so a SQLite database file must be provided.  By default, "database.sqlite" \
is used; an alternative database file can be specified using the -d option.')
argp.add_argument('-d', '--dbconf', help='the SQLite database file ("database.sqlite" by default)')
argp.add_argument('-l', '--resolvers', help=citehelpstr)
argp.add_argument('-c', '--comptaxoid', type=int, help='the ID of a taxonomy to check for name citation \
data (-1 [=none] by default)')
argp.add_argument('infile', help='the CSV taxonomy configuration file')
argp.set_defaults(dbconf='database.sqlite', resolvers='all', comptaxoid=-1)
args = argp.parse_args()

# Get a cursor for the taxonomy database.
try:
    pgcur = taxodatabase.getDBCursor(args.dbconf)
except ConfigError as e:
    exit('\n' + str(e) + '\n')

# Attempt to read the configuration file and parse the input taxonomy CSV file.
taxoconfig = TaxonomyConfig()
taxoparser = CSVTaxonomyParser()
try:
    taxoconfig.read(args.infile)
    print 'Parsing input CSV taxonomy file...'
    taxonomyroot = taxoparser.parseCSV(taxoconfig, pgcur)
    print 'done.'
except (ConfigError, TaxoCSVError) as e:
    exit('\n' + str(e) + '\n')

# Get the citation configuration information and create the Citation object.
citation = Citation(*taxoconfig.getCitationSettings())


# Get the taxonomy information from the configuration file.
taxonomyid, taxonomyname, ismaster = taxoconfig.getTaxonomySettings()

# Create the Taxonomy object.
taxonomy = Taxonomy(taxonomyid, taxonomyname, ismaster, citation, taxonomyroot)
print(taxonomyid, taxonomyname, taxonomyroot)
# Make sure that the Taxonomy is linked to the MOL backbone taxonomy.
#print '\nLinking taxonomy to the MOL backbone taxonomy...'
#if not(taxonomy.linkToBackbone(pgcur)):
#    exit('\nError:\n  Unable to link the taxonomy to the MOL backbone taxonomy.\n')
#print 'done.'
#taxonomy.getBackboneTaxonomy().printAll()

useres = args.resolvers
if useres != 'none':
    if useres != 'all':
        try:
            resindex = int(useres)
            resolvers = [resolvers[resindex]]
        except (ValueError, IndexError):
            exit('\nError: "' + useres + '" is an invalid name citations resolver index.\n' + 
                    'To use a specific resolver, you must provide an integer from 0 to '
                    + str(len(resolvers) - 1) + '.\n')

    # Try to retrieve citation information for the names in the taxonomy.
    # Note that we could include the backbone taxonomy in the name resolution process,
    # but we don't currently have a good way to get citation data for high-level taxa
    # anyway, so we just focus on the target taxonomy.
    for resolver in resolvers:
        print '\nBeginning name citation resolution with ' + resolver.getSourceDescription() + '...'
        if args.comptaxoid > -1:
            resolver.setComparisonTaxonomy(args.comptaxoid)
            resolver.setSkipIfExisting(True)
        resolver.resolve(pgcur, taxonomy.roottaxon)
        print 'name resolution with ' + resolver.getSourceDescription() + ' finished.'

# Make sure the taxonomy is persisted to the database.
print '\nPersisting taxonomy to the database...'
taxonomy.persist(pgcur, True)
print 'finished.'

totalrows, totaltaxa = taxoparser.getStats()
print '\nProcessed', totalrows, 'CSV file rows containing', totaltaxa, 'unique taxa.\n'

