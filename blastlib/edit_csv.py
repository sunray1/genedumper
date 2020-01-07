#!/usr/bin/env python
# Runs step 1
# Converts .xml to .csv and contains some filtering based on evalue and query coverage
# Pulls out species names from blast definitions as blast does not return the species individually
import csv, sys
import xml.etree.ElementTree as ET

class edit():

	def __init__(self, file):
		self.file = file
		baselist = self.file.rsplit(".", 1)
		basename = baselist[0]
		self.basename = basename

	def xmltocsv(self, taxdb):
		root = ET.parse(self.file).getroot()
		#  Get the length of the query sequence.
		querylen = int(root.find('BlastOutput_query-len').text)
		fieldnames = ['Name_num', 'Decision', 'num', 'Gene_name', 'GI', 'def', 'accession', 'hit_length',
					  'bit_score', 'evalue', 'q_from', 'q_to', 'h_from', 'h_to', 'identity',
					  'align_length', 'qseq', 'hseq', 'Species', 'genus', 'epithet',
					  'mismatches', 'qcov', 'tc_id']
		fout = open(self.basename + ".csv", 'w')
		writer = csv.DictWriter(fout, fieldnames)
		writer.writeheader()
		# Process each database hit.
		hits = root.find('BlastOutput_iterations/Iteration/Iteration_hits')
		row = {}
		for hit in hits:
			# Get the sequence information.
			idstr = hit.find('Hit_id').text
			idelems = idstr.split('|')
			GI = idelems[1]
			accession = idelems[3]
			hitlen = int(hit.find('Hit_len').text)
			defstr = hit.find('Hit_def').text
			num = int(hit.find('Hit_num').text)
			defstr.replace("\'", "\'\'")
				
			# Get the taxon information.
			pull_name_out = pull_names(defstr, accession, GI, taxdb)
			if len(pull_name_out) == 4:
				Species, genus, epithet, decision = pull_name_out
			elif len(pull_name_out) == 5:
				Species, genus, epithet, decision, GI = pull_name_out				
			# Get the first Hsp element.
			hsp = hit.find('Hit_hsps/Hsp[1]')
				
			# Get the alignment details.
			score = int(hsp.find('Hsp_score').text)
			bitscore = float(hsp.find('Hsp_bit-score').text)
			evalue = float(hsp.find('Hsp_evalue').text)
			alignlen = int(hsp.find('Hsp_align-len').text)
			identity = int(hsp.find('Hsp_identity').text)
			gaps = int(hsp.find('Hsp_gaps').text)
			q_from = int(hsp.find('Hsp_query-from').text)
			q_to = int(hsp.find('Hsp_query-to').text)
			h_from = int(hsp.find('Hsp_hit-from').text)
			h_to = int(hsp.find('Hsp_hit-to').text)
			qseq = hsp.find('Hsp_qseq').text
			hseq = hsp.find('Hsp_hseq').text
			if q_to > q_from:
				mismatch = q_to - q_from + 1 - identity
				qcov = float(q_to - q_from + 1) / querylen
			elif q_to < q_from:
				mismatch = q_from - q_to + 1 - identity
				qcov = float(q_from - q_to + 1) / querylen
			
			# Prepare the CSV output row and write it to the file,
			# provided it meets the acceptance criteria.
			row['num'] = num
			row['Gene_name'] = self.basename
			row['GI'] = GI
			row['def'] = defstr
			row['accession'] = accession
			row['hit_length'] = hitlen
			row['bit_score'] = bitscore
			row['evalue'] = evalue
			row['q_from'] = q_from
			row['q_to'] = q_to
			row['h_from'] = h_from
			row['h_to'] = h_to
			row['identity'] = identity
			row['align_length'] = alignlen
			row['qseq'] = qseq
			row['hseq'] = hseq
			row['Species'] = Species
			row['genus'] = genus
			row['epithet'] = epithet
			row['qcov'] = qcov
			row['mismatches'] = mismatch
			row['tc_id'] = 0
			row['Name_num'] = self.basename + "_" + str(num)
			row['Decision'] = decision
			if evalue < 0.001 and decision != "delete":
				writer.writerow(row)
				
#referenced in above function
def pull_names(defstr, accession, origGI, taxdb):
	import sqlite3
	conn = sqlite3.connect(taxdb)
	c = conn.cursor()
	#defstr = "Gynostemma 'burmanicum var. molle' isolate PT5203MT04 internal transcribed spacer 1, partial sequence; 5.8S ribosomal RNA gene, complete sequence; and internal transcribed spacer 2, partial sequence >gi|523928652|gb|KF269117.1| Gynostemma 'burmanicum var. molle' isolate PT5203MT05 internal transcribed spacer 1, partial sequence; 5.8S ribosomal RNA gene, complete sequence; and internal transcribed spacer 2, partial sequence"
	defstr = ">gi||||"+defstr
	#print(defstr)
	GI = ""
	matchbool = []
	decision = ""
	species ="Unknown species"
	#accession = "ehs"
	removelist = ["cf.", "nr.", " aff. ", " hybrid ", " form ", " n. sp.", " sp. n.", " Cf. ", " x ", " X ", " sp. ", " gen. "]
	for i in defstr.split(">gi")[1:]:
		#loops through till finds sp in taxonomy
		if species == "Unknown species":
			name_id = ""
			defstr_sub = i.split("|")[4]
			defin = defstr_sub.split()
			#print("1")
			if any(x in defstr_sub for x in removelist):
				species = "Ambiguous species"
				decision = "Ambiguous species/not chosen"
				#print("2")
				break
			else:
				#print("3")
				species = defin[0] + " " + defin[1]
				species = species.replace("\"", "")
				species = species.replace("'", "")
				species = species.replace("(", "")
				species = species.replace(")", "")
				for iter in c.execute("SELECT name_id FROM names WHERE namestr = '"+species+"'"):
					#will only iterate through if species exists
					GI = i.split("|")[1]
					if GI == "":
						GI = origGI
					decision = ""
					#print("4")
				if GI != "":
					#print("5")
					break
				else:
					#print("6")
					decision = "Species not in taxonomy/not chosen"
		
	# if boolian = True:
	# sometimes definitions look like 'PREDICTED: Genus Species'
	if accession[0] == "X" and accession[2] == "_":
		decision = "delete"

	try:
		genus = species.split()[0].replace("'", "")
		epithet = species.split()[1].replace("'", "")
	except:
		species = "Unknown species"
		genus = species.split()[0].replace("'", "")
		epithet = species.split()[1].replace("'", "")
		decision = "Unknown species/not chosen"
	species = species.replace("'", "")
	#print(species, decision, GI)
	if GI != "":
		return(species, genus, epithet, decision, GI)
	else:
		return(species, genus, epithet, decision)



