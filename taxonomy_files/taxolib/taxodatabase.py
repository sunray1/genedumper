
#import psycopg2 as ppg2
import sqlite3
from ConfigParser import RawConfigParser
from taxoconfig import ConfigError


def _getPostgresDBCursor(conffile):
    # Read the database connection settings from the configuration file.
    cp = RawConfigParser()
    res = cp.read(conffile)
    if len(res) == 0:
        raise ConfigError('The database configuration file ' + conffile + ' could not be read.')

    dbhost = cp.get('database', 'host')
    dbport = cp.getint('database', 'port')
    dbuser = cp.get('database', 'user')
    dbpw = cp.get('database', 'pw')
    dbname = cp.get('database', 'dbname')
    dbschema = cp.get('database', 'schema')

    # Set up the database connection, get a connection cursor, and set the default schema.
    conn = ppg2.connect(host=dbhost, port=dbport, user=dbuser, password=dbpw, database=dbname)
    pgcur = conn.cursor()
    pgcur.execute('SET search_path TO ' + dbschema)

    return pgcur

def _getSQLiteDBCursor(dbfile):
    conn = sqlite3.connect(dbfile)
    slcur = conn.cursor()

    return slcur

def getDBCursor(filein):
    return _getSQLiteDBCursor(filein)

