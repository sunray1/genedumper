#!/usr/bin/env python
# This script is designed to run step 6 - searching for misspellings in Species names. Unfortunately
# people don't misspell the same way every time so some of this must be done by hand. This script will
# produce a list of changed Species and suggested changes called 'spell.txt'
def spelling(taxdb, blastdb):

	import sqlite3
	name_num_dic = {}
	conn = sqlite3.connect(taxdb)
	c = conn.cursor()
	c.execute("ATTACH '" + blastdb + "' as 'db'")
	for iter in c.execute("SELECT Name_num, Species FROM blast WHERE tc_id IS NULL and epithet != 'sp.' and epithet != 'gen.' and epithet NOT LIKE '%/%';"):
		if str(iter[1]) in name_num_dic.keys():
			list1 = list(name_num_dic[str(iter[1])])
			list1.append(str(iter[0]))
			name_num_dic[str(iter[1])] = list1
		else:
			name_num_dic[str(iter[1])] = [str(iter[0])]
	with open('spell.txt', 'w') as o:
		for key in name_num_dic.keys():
			update = False
			for num in name_num_dic[key]:
				if num == name_num_dic[key][0]:
					#First, lets check to see if species are misnamed by using an 'a' instead of 'us' or vice versa
					if key[-2:] == 'us':
						organism = []
						name_ch = key[:-2] + 'a'
						for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt, taxon_concepts tc, ranks r WHERE n.namestr LIKE '" + name_ch + "%' and n.name_id = ntt.name_id and ntt.tc_id = tc.tc_id and tc.rank_id=r.rank_id and r.namestr = 'Species'"):
							organism.append(iter[0])
							organism = set(organism)
							organism = list(organism)
						if len(organism) == 1:
							organism = organism[0]
							c.execute("UPDATE blast SET Species='" + organism + "' WHERE Name_num='" + num + "';")
							c.execute("UPDATE blast SET genus='" + organism.split()[0] + "' WHERE Name_num='" + num + "';")
							c.execute("UPDATE blast SET epithet='" + organism.split()[1] + "' WHERE Name_num='" + num + "';")
							o.write(key + '\t' + str(name_num_dic[key]).replace("'", "").strip("[").strip("]") + '\t' + organism + '\tchanged\n')
							update = True
					elif key[-1:] == 'a':
						organism = []
						name_ch = key[:-1] + 'us'
						for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt, taxon_concepts tc, ranks r WHERE n.namestr LIKE '" + name_ch + "%' and n.name_id = ntt.name_id and ntt.tc_id = tc.tc_id and tc.rank_id=r.rank_id and r.namestr = 'Species'"):
							organism.append(iter[0])
							organism = set(organism)
							organism = list(organism)
						if len(organism) == 1:
							organism = organism[0]
							c.execute("UPDATE blast SET Species='" + organism + "' WHERE Name_num='" + num + "';")
							c.execute("UPDATE blast SET genus='" + organism.split()[0] + "' WHERE Name_num='" + num + "';")
							c.execute("UPDATE blast SET epithet='" + organism.split()[1] + "' WHERE Name_num='" + num + "';")
							o.write(key + '\t' + str(name_num_dic[key]).replace("'", "").strip("[").strip("]") + '\t' + organism + '\tchanged\n')
							update = True
					# Now lets check to see if species are misnamed due to an 'i' instead of an 'ii'
					elif key[-1:] == 'i':
						organism = []
						name_ch = key[:-1] + 'ii'
						for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt, taxon_concepts tc, ranks r WHERE n.namestr LIKE '" + name_ch + "%' and n.name_id = ntt.name_id and ntt.tc_id = tc.tc_id and tc.rank_id=r.rank_id and r.namestr = 'Species'"):
							organism.append(iter[0])
							organism = set(organism)
							organism = list(organism)
						if len(organism) == 1:
							organism = organism[0]
							c.execute("UPDATE blast SET Species='" + organism + "' WHERE Name_num='" + num + "';")
							c.execute("UPDATE blast SET genus='" + organism.split()[0] + "' WHERE Name_num='" + num + "';")
							c.execute("UPDATE blast SET epithet='" + organism.split()[1] + "' WHERE Name_num='" + num + "';")
							o.write(key + '\t' + str(name_num_dic[key]).replace("'", "").strip("[").strip("]") + '\t' + organism + '\tchanged\n')
							update = True
					else:
						# Here we just check for different endings in general by checking to see if the beginnings of any of the names match anything in the
						# taxonomy database. If there are multiple hits, it won't be changed but it will be noted in 'spell.txt' and can be changed by hand
						organism = []
						name_ch = key[:-2]
						for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt, taxon_concepts tc, ranks r WHERE n.namestr LIKE '" + name_ch + "%' and n.name_id = ntt.name_id and ntt.tc_id = tc.tc_id and tc.rank_id=r.rank_id and r.namestr = 'Species'"):
							organism.append(iter[0])
							organism = set(organism)
							organism = list(organism)
						if len(organism) == 1:
							organism = organism[0]
							c.execute("UPDATE blast SET Species='" + organism + "' WHERE Name_num='" + num + "';")
							c.execute("UPDATE blast SET genus='" + organism.split()[0] + "' WHERE Name_num='" + num + "';")
							c.execute("UPDATE blast SET epithet='" + organism.split()[1] + "' WHERE Name_num='" + num + "';")
							o.write(key + '\t' + str(name_num_dic[key]).replace("'", "").strip("[").strip("]") + '\t' + organism + '\tchanged\n')
							update = True
						elif len(organism) > 1:	
							o.write(key + '\t' + str(name_num_dic[key]).replace("'", "").strip("[").strip("]") + '\t' + str(organism) + '\tmultiple hits\n')
					spes = key.split()[1]
					#just change these by hand
					# This tries to deal with mispellings in general by searching for species that have names like the beginnings and ends of null species.
					# However, because the results from this are so variable, they are noted in 'spell.txt' and must be judged and changed by hand
					if len(spes) > 4 and update != True:
						organism = []
						name_ch = spes[:2] + '%' + spes[-2:]
						for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt, taxon_concepts tc, ranks r WHERE n.namestr LIKE '" + key.split()[0] + " " + name_ch + "' and n.name_id = ntt.name_id and ntt.tc_id = tc.tc_id and tc.rank_id=r.rank_id and r.namestr = 'Species'"):
							organism.append(iter[0])
							organism = set(organism)
							organism = list(organism)
						if len(organism) == 1:
							organism = organism[0]
							o.write(key + '\t' + str(name_num_dic[key]).replace("'", "").strip("[").strip("]") + '\t' + str(organism) + '\thand\n')
					else:
						pass
				elif num != name_num_dic[key][0] and update == True:
					c.execute("UPDATE blast SET Species='" + organism + "' WHERE Name_num='" + num + "';")
					c.execute("UPDATE blast SET genus='" + organism.split()[0] + "' WHERE Name_num='" + num + "';")
					c.execute("UPDATE blast SET epithet='" + organism.split()[1] + "' WHERE Name_num='" + num + "';")
	conn.commit()
	conn.close()
