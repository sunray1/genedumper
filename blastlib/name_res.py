#!/usr/bin/env python
# This script is used to update the tc_ids once the Species names are changed. This can be also be used after the Species
# are changed by hand
import sys, sqlite3

def find_id(taxonomy, blastfile):	
	conn = sqlite3.connect(blastfile)
	c = conn.cursor()
	c.execute("ATTACH '" + taxonomy + "' as 'tax'")
	try:
		c.execute("ALTER TABLE blast ADD tc_id int")
	except:
		pass
	c.execute("UPDATE blast SET tc_id = (SELECT tc.tc_id FROM taxon_concepts tc, names_to_taxonconcepts ntt, names n WHERE tc.tc_id = ntt.tc_id AND ntt.name_id = n.name_id AND n.namestr = blast.Species);")
	conn.commit()
	conn.close()