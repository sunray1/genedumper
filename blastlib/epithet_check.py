#!/usr/bin/env python
# Runs step 5
# This script searches the names of all the tc_id=null species against the taxonomy database, looking
# for possible valid names- it does this in two ways
# First the Genus and epithet are used to search for matching subspecies with matching Genera, but
# the epithet is instead the subspecific epithet. If there is a hit, the Species for all the blast
# hits containing the previous name are changed.

# Second the epithets are used to search for all Species with the same epithet. If this returned list
# contains species with the same tc_id and there is only one valid name in the list, the Species for
# all the blast hits containing the previous name are changed.

# Changes are noted in 'epithet.txt'
def epithet(taxdb, blastdb):

	import sqlite3
	name_num_dic = {}
	conn = sqlite3.connect(taxdb)
	c = conn.cursor()
	conn.text_factory = str
	c.execute("ATTACH '" + blastdb + "' as 'db'")
	for iter in c.execute("SELECT Name_num, epithet, genus FROM blast WHERE tc_id IS NULL and epithet != 'sp.' and epithet != 'gen.' and epithet NOT LIKE '%/%';"):
		taxname = str(iter[2]) + " " + str(iter[1])
		if taxname in name_num_dic.keys():
			list1 = list(name_num_dic[taxname])
			list1.append(str(iter[0]))
			name_num_dic[taxname] = list1
		else:
			name_num_dic[taxname] = [str(iter[0])]
	with open('epithet.txt', 'w') as o:
		for key in name_num_dic.keys():
			update = False
			for num in name_num_dic[key]:
				if num == name_num_dic[key][0]:
					name = key.split()[1]
					genus = key.split()[0]
					species = []
					tc_id = []
					validity = []
					organism = []
#					Searches subspecies
#					there are no intances where both conditions are true so won't overlap(at least for BN)
					for iter in c.execute("SELECT tc.tc_id, n.namestr FROM names n, names_to_taxonconcepts ntt, taxon_concepts tc, ranks r WHERE n.namestr LIKE '" + genus + " % " + name + "' and n.name_id = ntt.name_id and ntt.tc_id = tc.tc_id and tc.rank_id=r.rank_id and r.namestr = 'Subspecies' and ntt.validity = 'valid';"):
						organism.append(iter[1])
						
					if len(organism) == 1:
						organism = organism[0]
						c.execute("UPDATE blast SET Species='" + organism + "' WHERE Name_num='" + num + "';")
						c.execute("UPDATE blast SET genus='" + organism.split()[0] + "' WHERE Name_num='" + num + "';")
						c.execute("UPDATE blast SET epithet='" + organism.split()[1] + "' WHERE Name_num='" + num + "';")
						o.write('Subspecies\t' + genus + ' ' + name +'\t' + str(name_num_dic[key]).replace("'", "").strip("[").strip("]") + '\t' + organism + '\n')
						update = True
#					Searches species
					else:
						for iter in c.execute("SELECT tc.tc_id, n.namestr, ntt.validity FROM names n, names_to_taxonconcepts ntt, taxon_concepts tc, ranks r WHERE n.namestr LIKE '% " + name + "' and n.name_id = ntt.name_id and ntt.tc_id = tc.tc_id and tc.rank_id=r.rank_id and r.namestr = 'Species';"):
							tc_id.append(iter[0])
							species.append(iter[1])
							validity.append(iter[2])
#							have to make set because sometimes there are duplicates in valids and syns - need to talk to Vijay about this
						tc_id = set(tc_id)
						if len(tc_id) == 1 and validity.count('valid') == 1:
							organism = species[validity.index('valid')]
							c.execute("UPDATE blast SET Species='" + organism + "' WHERE Name_num='" + num + "';")
							c.execute("UPDATE blast SET genus='" + organism.split()[0] + "' WHERE Name_num='" + num + "';")
							c.execute("UPDATE blast SET epithet='" + organism.split()[1] + "' WHERE Name_num='" + num + "';")
							o.write('Species\t' + genus + ' ' + name + '\t' + str(name_num_dic[key]).replace("'", "").strip("[").strip("]") + '\t' + organism + '\n')
							update = True
				elif num != name_num_dic[key][0] and update == True:
					c.execute("UPDATE blast SET Species='" + organism + "' WHERE Name_num='" + num + "';")
					c.execute("UPDATE blast SET genus='" + organism.split()[0] + "' WHERE Name_num='" + num + "';")
					c.execute("UPDATE blast SET epithet='" + organism.split()[1] + "' WHERE Name_num='" + num + "';")
						
	conn.commit()
	conn.close()