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
	for iter in c.execute("SELECT tc_id , count(*) FROM (SELECT tc_id, Gene_name FROM blast WHERE tc_id NOT NULL GROUP BY tc_id, Gene_name) GROUP BY tc_id;"):
		iddic[int(iter[0])] = int(iter[1])
		
	#get % with 0 markers
	
	for iter in c.execute("SELECT tc_id , count(*) FROM (SELECT tc_id, Gene_name FROM blast WHERE tc_id NOT NULL GROUP BY tc_id, Gene_name) GROUP BY tc_id;"):
		iddic[int(iter[0])] = int(iter[1])

	gene = []
	gene_amount = []
	#get top 3 genes that have tc_ids
	for iter in c.execute("SELECT Gene_name, count(*) as count from (SELECT Gene_name, tc_id from blast WHERE tc_id NOT NULL GROUP BY Gene_name, tc_id) GROUP BY Gene_name ORDER BY count DESC;"):
		gene.append(iter[0])
		gene_amount.append(int(iter[1]))

	loci_num = sorted(zip(gene_amount, gene), reverse=True)
	
	decision = []
	decision_amount = []
	#get top 3 genes that have tc_ids
	for iter in c.execute("SELECT decision, count(*) as count FROM blast WHERE tc_id not null GROUP BY decision ORDER BY count DESC;"):
		decision.append(iter[0])
		decision_amount.append(int(iter[1]))
	
	decision_num = sorted(zip(decision_amount, decision), reverse=True)

	twentyfiveperc = []
	fiftyperc = []
	seventyfiveperc = []
	ninetyfiveperc = []
	for i, k in enumerate(gene_amount):
		gene_perc = float(k/total_num_species)
		if gene_perc > .25:
			twentyfiveperc.append(gene[i])
		if gene_perc >= .5:
			fiftyperc.append(gene[i])
		if gene_perc >= .75:
			seventyfiveperc.append(gene[i])
		if gene_perc >= .95:
			ninetyfiveperc.append(gene[i])

	
	
	perc_spec_markers_zero = int(round(float(total_num_species-len(iddic))/float(total_num_species), 2)*100)
	## get tree data
	rank_count = 0
	for iter in c.execute("select rank_id, count(*) as count from taxon_concepts group by rank_id order by rank_id DESC"):
			rank_id = iter[0]
			rank_count = int(iter[1])
			if rank_count == 1:
				break

	rangelist = list(range(rank_id, 135, 5))
	rangelist = str(list(map(lambda x:str(x), rangelist))).replace("[", "(").replace("]", ")")
	namestr = []
	parent_id_list = []
	tc_id_list = []
	for iter in c.execute("SELECT names.namestr, taxon_concepts.parent_id, taxon_concepts.tc_id FROM\
						  names JOIN names_to_taxonconcepts on names_to_taxonconcepts.name_id=names.name_id\
						  JOIN taxon_concepts on taxon_concepts.tc_id = names_to_taxonconcepts.tc_id\
						  WHERE taxon_concepts.rank_id in "+rangelist+" and names_to_taxonconcepts.validity = 'valid'"):
		namestr.append(iter[0])
		parent_id_list.append(iter[1])
		tc_id_list.append(iter[2])
	parent_id_list_str = str(parent_id_list).replace("[", "(").replace("]", ")")
	parent_namestr_dic = {}
	parent_namestr = []
	for iter in c.execute("SELECT names.namestr, taxon_concepts.tc_id FROM names JOIN names_to_taxonconcepts\
						  on names_to_taxonconcepts.name_id=names.name_id JOIN taxon_concepts\
						  on taxon_concepts.tc_id = names_to_taxonconcepts.tc_id WHERE\
						  taxon_concepts.tc_id IN " + parent_id_list_str + " and names_to_taxonconcepts.validity = 'valid'"):
		parent_namestr_dic[int(iter[1])] = iter[0]
	#we are going to remove the first one since it is the base taxon anyway
	parent_namestr_dic[0] = 'null'
	for i in parent_id_list:
		parent_namestr.append(parent_namestr_dic[i])
	genecountdic = {}
	for iter in c.execute("SELECT tc_id, count(tc_id) from (SELECT tc_id, Gene_name from blast WHERE tc_id NOT NULL GROUP BY Gene_name, tc_id) group by tc_id"):
		genecountdic[iter[0]] = int(iter[1])
	gene_counts = []
	for i in tc_id_list:
		try:
			gene_counts.append(genecountdic[i])
		except KeyError:
			gene_counts.append(0)
	data = list(zip(namestr, parent_namestr, [1] * len(namestr), gene_counts))

	#get gene_counts using SELECT tc_id, count(tc_id) from blast WHERE tc_id NOT NULL GROUP BY Gene_name, tc_id
	#get names from tc_id = SELECT namestr FROM names JOIN names_to_taxonconcepts ON names.name_id=names_to_taxonconcepts.name_id WHERE names_to_taxonconcepts.tc_id = '120' AND names_to_taxonconcepts.validity = 'valid'
	
			
	with open("info.html", "w") as o:
		o.write("<!DOCTYPE html>\n")
		o.write("<HTML><HEAD>\n")
		o.write("<TITLE>Hero Stats</TITLE>\n")
		
		o.write("<style>\n\
* {\n\
  box-sizing: border-box;\n\
}\n\
\
.box {\n\
  float: left;\n\
  width: 50%;\n\
  padding: 50px;\n\
}\n\
\
.clearfix::after {\n\
  content: \"\";\n\
  clear: both;\n\
  display: table;\n\
}\n\
</style>\n")
		o.write("</HEAD>\n")
		o.write("<BODY>\n")
		o.write("<h1>Hero Statistics</h1>\n")
		o.write("<p>Total number of species: " + str(total_num_species) + "<br>\n")
		o.write("Total number of species with at least one marker: " + str(tc_id_inblast) + " (" + perc_species + "%)<br>\n")
		o.write("<br>")
		o.write("Number of loci with > 25% species: " + str(len(twentyfiveperc)) + "<br>\n")
		o.write("Number of loci with > 50% species: " + str(len(fiftyperc)) + "<br>\n")
		o.write("Number of loci with > 75% species: " + str(len(seventyfiveperc)) + "<br>\n")
		o.write("Number of loci with > 95% species: " + str(len(ninetyfiveperc)) + "<br>\n")
		o.write("</p>\n")
		
		o.write("\n")
		o.write("\n")
		o.write("\n")
		o.write("<div id=\"piechart\"></div>\n")
		o.write("<div class=\"clearfix\">\n")
		o.write("<div class=\"box\" style=\"\">\n")
		o.write("<div id=\"barchart_values\" style=\"\"></div>\n")
		o.write("</div>\n")
		o.write("<div class=\"box\" style=\"\">\n")
		o.write("<div id=\"columnchart_values\" style=\"\"></div>\n")
		o.write("</div>\n")
		o.write("</div>\n")
		o.write("<div id=\"chart_div\" style=\"height: 500px;\"></div>\n")
		o.write("</BODY>\n")
		
		

		
		#piechart
		o.write("<script type=\"text/javascript\" src=\"https://www.gstatic.com/charts/loader.js\"></script>\n")
		o.write("<script type=\"text/javascript\"> \n // Load google charts \n google.charts.load('current', {'packages':['corechart']});\n google.charts.setOnLoadCallback(drawChart);\n")
		o.write("function drawChart() {\
			var data = google.visualization.arrayToDataTable([\n\
			['Number of Markers', 'Number of species with marker'],\n")
		o.write("['0', "+str(total_num_species-len(iddic))+"],\n")
		for i in range(num_markers):
			perc_spec_markers = list(iddic.values()).count(i+1)
			o.write("['" + str(i+1) + "', " + str(perc_spec_markers) + "],\n")
		o.write("]);\n")
		o.write("var options = {'title':'Percent of species with x marker(s)', 'width':550, 'height':400};\n")
		o.write("var chart = new google.visualization.PieChart(document.getElementById('piechart'));\n")
		o.write("chart.draw(data, options);}\n</script>\n")
		o.write("\n\n\n")
		
		#loci bar chart
		o.write("<script type=\"text/javascript\" src=\"https://www.gstatic.com/charts/loader.js\"></script>\n")
		o.write("<script type=\"text/javascript\"> \n // Load google charts \n google.charts.load('current', {'packages':['corechart']});\n google.charts.setOnLoadCallback(drawChart);\n")
		o.write("function drawChart() {\
			var data = google.visualization.arrayToDataTable([\n\
			['Locus', 'Number of species', { role: \"style\" } ],\n")
		for i, k in loci_num: 
			o.write("['" + str(k) + "', " + str(i) + ", \"dark blue\"],\n")
		o.write("]);\n")
			
		o.write("var view = new google.visualization.DataView(data);\n\
					view.setColumns([0, 1,\n\
                       { calc: \"stringify\",\n\
                         sourceColumn: 1,\n\
                         type: \"string\",\n\
                         role: \"annotation\" },\n\
                       2]);\n\
				")
		o.write("var options = {\
        title: \"Number of species per locus\",\
        width: 600,\
        height: 400,\
        bar: {groupWidth: \"75%\"},\
        legend: { position: \"none\" },\
      };")
		
		o.write("var chart = new google.visualization.BarChart(document.getElementById(\"barchart_values\"));\n\
			chart.draw(view, options);\n\
			}\
			</script>\n")
		o.write("\n\n\n")
		###decision chart

		o.write("<script type=\"text/javascript\" src=\"https://www.gstatic.com/charts/loader.js\"></script>\n")
		o.write("<script type=\"text/javascript\"> \n // Load google charts \n\tgoogle.charts.load('current', {'packages':['corechart']});\n\tgoogle.charts.setOnLoadCallback(drawChart);\n")
		o.write("\tfunction drawChart() {\n\
			\tvar data = google.visualization.arrayToDataTable([\n\
			\t['Decision', 'Number of hits', { role: \"style\" } ],\n")
		for i, k in decision_num: 
			o.write("['" + str(k) + "', " + str(i) + ", \"dark red\"],\n")
		o.write("]);\n")
			
		o.write("      var view = new google.visualization.DataView(data);\n\
						view.setColumns([0, 1,\n\
                       { calc: \"stringify\",\n\
                         sourceColumn: 1,\n\
                         type: \"string\",\n\
                         role: \"annotation\" },\n\
                       2]);\n")
		o.write("      var options = {\n\
        title: \"Decisions amounts\",\n\
        width: 600,\n\
        height: 400,\n\
        bar: {groupWidth: \"75%\"},\n\
        legend: { position: \"none\" },\n\
      };\n")
		
		o.write("var chart = new google.visualization.ColumnChart(document.getElementById(\"columnchart_values\"));\n\
			chart.draw(view, options);\n\
			}\n\
			</script>\n")
		o.write("\n\n\n")
		
		##tree chart
		
		o.write("<script type=\"text/javascript\" src=\"https://www.gstatic.com/charts/loader.js\"></script>\n")
		o.write("<script type=\"text/javascript\">\n\
      google.charts.load('current', {'packages':['treemap']});\n\
      google.charts.setOnLoadCallback(drawChart);\n")
		o.write("function drawChart() {\n\
        var data = google.visualization.arrayToDataTable([\n")
		o.write("\t['Taxon', 'Parent taxon', 'size', 'Number of loci (color)'],\n")
		o.write("\t['" + data[0][0] + "', null, 1, 0],\n")
		for i in data[1:]:
			o.write("\t" + str(i).replace("(", "[").replace(")", "]") + ",\n")
		###data goes here
		o.write("\t]);\n")
		o.write("\n")
		o.write("tree = new google.visualization.TreeMap(document.getElementById('chart_div'));\n")
		o.write("tree.draw(data, {\n\
          minColor: '#ff0000',\n\
          midColor: '#ffff00',\n\
          maxColor: '#00ff00',\n\
          headerHeight: 15,\n\
          fontColor: 'black',\n\
          showScale: true,\n\
          generateTooltip: showFullTooltip\n\
        });\n\n")
		o.write("  function showFullTooltip(row, size) {\n\
			return '<div style=\"background:#fd9; padding:10px; border-style:solid\">' +\n\
           '<span style=\"font-family:Courier\"><b>' + data.getValue(row, 0) +\n\
           '</b>, ' + data.getValue(row, 1) + ', ' + data.getValue(row, 2) +\n\
           ', ' + data.getValue(row, 3) + '</span><br>' +\n\
           'Value of this cell: ' + data.getValue(row, 2) + '<br>' + \n\
           data.getColumnLabel(2) +\n\
           ' (total value of this cell and its children): ' + size + '<br>' +\n\
			data.getColumnLabel(3) + ': ' + data.getValue(row, 3) + ' </div>';\n\
			}\n\n")
		o.write ("}\n")
		o.write("</script>\n")
		o.write("</HTML>\n")

	conn.commit()
	conn.close()
			
