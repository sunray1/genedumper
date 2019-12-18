#!/usr/bin/env python
# This script is designed to create an HTML file containing information about the BLAST table
import sqlite3
def html(taxonomy, blastfile):
	conn = sqlite3.connect(taxonomy)
	c = conn.cursor()
	c.execute("ATTACH '" + blastfile + "' as 'blastf'")
	
	
	#paragraph 1
	for iter in c.execute("SELECT COUNT(*) FROM (SELECT tc_id FROM blast WHERE tc_id NOT NULL GROUP BY tc_id)"):
		tc_id_inblast = int(iter[0])
		
	for iter in c.execute("SELECT COUNT(*) FROM (SELECT n.namestr\
						   FROM names n JOIN names_to_taxonconcepts ntt ON ntt.name_id = n.name_id\
						   JOIN taxon_concepts tc ON tc.tc_id = ntt.tc_id\
						   JOIN ranks r ON r.rank_id = tc.rank_id WHERE r.namestr = 'Species' AND ntt.validity='valid' GROUP BY n.namestr)"):
		total_num_species = int(iter[0])
	perc_species = str(int(round(float(tc_id_inblast)/float(total_num_species), 2)*100))
	
	
	#paragraph 2
	
	for iter in c.execute("SELECT COUNT(*) FROM (SELECT Gene_name FROM blast GROUP BY Gene_name)"):
		num_markers = int(iter[0])

	iddic = {}
	#make dic with tc_ids and gene amounts in blast
	for iter in c.execute("SELECT tc_id , count(*) FROM (SELECT tc_id, Gene_name FROM blast WHERE tc_id NOT NULL GROUP BY tc_id, Gene_name) GROUP BY tc_id"):
		iddic[int(iter[0])] = int(iter[1])
		
	#get % with 0 markers
	
	perc_spec_markers_zero = int(round(float(total_num_species-len(iddic))/float(total_num_species), 2)*100)
			
	with open("info.html", "w") as o:
		o.write("<!DOCTYPE html>\n")
		o.write("<HTML><HEAD>\n")
		o.write("<TITLE>Hero Stats</TITLE></HEAD>\n")
		o.write("<BODY>\n")
		
		o.write("<p>Total number of species: " + str(total_num_species) + "<br>\n")
		o.write("Total number of species with at least one marker: " + str(tc_id_inblast) + " (" + perc_species + "%)\n")
		o.write("</p>")
		
		o.write("Percent of species with 0 marker(s): " + str(perc_spec_markers_zero) + "%<br>\n")
		#percent of species with i marker NOT percent of species with AT LEAST i markers
		for i in range(num_markers):
			perc_spec_markers = int(round(float(list(iddic.values()).count(i+1))/float(total_num_species), 2)*100)
			o.write("Percent of species with " + str(i+1) + " marker(s): " + str(perc_spec_markers) + "%<br>\n")
		
		o.write("<p>")
		o.write("</p>")
		o.write("</BODY>\n")
		o.write("</HTML>\n")
		


# 	for iter in c.execute("SELECT n.namestr, tc.tc_id FROM names n JOIN names_to_taxonconcepts ntt ON ntt.name_id = n.name_id JOIN taxon_concepts tc ON tc.tc_id = ntt.tc_id JOIN ranks r ON r.rank_id = tc.rank_id WHERE r.namestr = '" + level + "' AND ntt.validity='valid' GROUP BY n.namestr;"):
# 		ranks.append(iter[0])
# 		tc_id.append(iter[1])
# 	#make dic with tc_ids and gene amounts in blast
# 	for iter in c.execute("SELECT tc_id, Gene_name FROM blast WHERE tc_id NOT NULL GROUP BY tc_id, Gene_name"):
# 		if int(iter[0]) in iddic.keys():
# 			list1 = list(iddic[int(iter[0])])
# 			list1.append(str(iter[1]))
# 			iddic[int(iter[0])] = list1
# 		else:
# 			iddic[int(iter[0])] = [str(iter[1])]
# 	for i in tc_id:
# 		idnums = [i]
# 		typelist = ["Family"]
# 		#make empty dictionary with genes and their amounts
# 		for iter in c.execute("SELECT Gene_name FROM blast GROUP BY Gene_name;"):
# 			genedic[iter[0]] = 0
# 		#make this into  two lists so it skips species  so dic[id] = species and if not 'species'	
# 		#gets list of all tc_ids in lower taxonomy for each valid name in level
# 		for n in idnums:
# 			if typelist[idnums.index(n)] != 'Species':
# 				for iter in c.execute("SELECT tc.tc_id, r.namestr FROM taxon_concepts tc, ranks r WHERE tc.rank_id=r.rank_id AND tc.parent_id =" + str(n) + ";"):
# 					idnums.append(iter[0])
# 					typelist.append(str(iter[1]))
# 		#for list of tc_ids for each level, check to see if in blast
# 		for n in idnums:
# 			if n in iddic.keys():
# 				for n in iddic[n]:
# 					genedic[n] = 1
# 		with open(level + "_stats.csv", "a") as o:
# 			o.write("\n" + ranks[tc_id.index(i)] + ",")
# 			for i in sorted(genedic):
# 				o.write(str(genedic[i]) + ",")
# 			o.write(str(sum(genedic.values())))
# 	# If the level is Species
# #	else:
# #		#make empty dictionary with genes and their amounts
# #		for i in tc_id:
# #			for iter in c.execute("SELECT Gene_name FROM blast GROUP BY Gene_name;"):
# #				genedic[iter[0]] = 0
# #			if i in iddic.keys():
# #				for n in iddic[i]:
# #					genedic[n] = 1
# #			with open(level + "_stats.csv", "a") as o:
# #				o.write("\n" + ranks[tc_id.index(i)] + ",")
# #				for i in sorted(genedic):
# #					o.write(str(genedic[i]) + ",")
# #				o.write(str(sum(genedic.values())))        

	conn.commit()
	conn.close()
			
