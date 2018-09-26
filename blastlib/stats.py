#!/usr/bin/env python
# This script is designed to recieve statistics about the blast results - specifically a matrix with
# the genes as columns and desired level (ie Species or Genera) as rows. Intersects are marked '1' if
# that level has that gene or a '0' if it does not. 
import sqlite3
def get_stats(taxonomy, blastfile, level):
	conn = sqlite3.connect(taxonomy)
	c = conn.cursor()
	c.execute("ATTACH '" + blastfile + "' as 'blastf'")
	with open(level + "_stats.csv", "w") as o:
		for iter in c.execute("SELECT Gene_name FROM blast GROUP BY Gene_name;"):
			o.write("," + iter[0])
	tc_id= []
	ranks= []
	genedic = {}
	iddic = {}
#first get list of every valid name and tc_id in that level
#need to make not case-sensitive
	for iter in c.execute("SELECT n.namestr, tc.tc_id FROM names n JOIN names_to_taxonconcepts ntt ON ntt.name_id = n.name_id JOIN taxon_concepts tc ON tc.tc_id = ntt.tc_id JOIN ranks r ON r.rank_id = tc.rank_id WHERE r.namestr = '" + level + "' AND ntt.validity='valid' GROUP BY n.namestr;"):
		ranks.append(iter[0])
		tc_id.append(iter[1])
	#make dic with tc_ids and gene amounts in blast
	for iter in c.execute("SELECT tc_id, Gene_name FROM blast WHERE tc_id NOT NULL GROUP BY tc_id, Gene_name"):
		if int(iter[0]) in iddic.keys():
			list1 = list(iddic[int(iter[0])])
			list1.append(str(iter[1]))
			iddic[int(iter[0])] = list1
		else:
			iddic[int(iter[0])] = [str(iter[1])]
	for i in tc_id:
		idnums = [i]
		typelist = ["Family"]
		#make empty dictionary with genes and their amounts
		for iter in c.execute("SELECT Gene_name FROM blast GROUP BY Gene_name;"):
			genedic[iter[0]] = 0
		#make this into  two lists so it skips species  so dic[id] = species and if not 'species'	
		#gets list of all tc_ids in lower taxonomy for each valid name in level
		for n in idnums:
			if typelist[idnums.index(n)] != 'Species':
				for iter in c.execute("SELECT tc.tc_id, r.namestr FROM taxon_concepts tc, ranks r WHERE tc.rank_id=r.rank_id AND tc.parent_id =" + str(n) + ";"):
					idnums.append(iter[0])
					typelist.append(str(iter[1]))
		#for list of tc_ids for each level, check to see if in blast
		for n in idnums:
			if n in iddic.keys():
				for n in iddic[n]:
					genedic[n] = 1
		with open(level + "_stats.csv", "a") as o:
			o.write("\n" + ranks[tc_id.index(i)] + ",")
			for i in sorted(genedic):
				o.write(str(genedic[i]) + ",")
			o.write(str(sum(genedic.values())))
	# If the level is Species
#	else:
#		#make empty dictionary with genes and their amounts
#		for i in tc_id:
#			for iter in c.execute("SELECT Gene_name FROM blast GROUP BY Gene_name;"):
#				genedic[iter[0]] = 0
#			if i in iddic.keys():
#				for n in iddic[i]:
#					genedic[n] = 1
#			with open(level + "_stats.csv", "a") as o:
#				o.write("\n" + ranks[tc_id.index(i)] + ",")
#				for i in sorted(genedic):
#					o.write(str(genedic[i]) + ",")
#				o.write(str(sum(genedic.values())))        

	conn.commit()
	conn.close()
			
