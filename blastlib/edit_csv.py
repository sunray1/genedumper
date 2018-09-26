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

	def xmltocsv(self):
		root = ET.parse(self.file).getroot()
		#  Get the length of the query sequence.
		querylen = int(root.find('BlastOutput_query-len').text)
		fieldnames = ['Name_num', 'num', 'Gene_name', 'GI', 'def', 'accession', 'hit_length',
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
			Species, genus, epithet = pull_names(defstr)
				
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
			if evalue < 0.001:
				writer.writerow(row)
				
#referenced in above function
def pull_names(defstr):
	# a lot of this is not generalized - need to figure out a way to change this
	removelist = ["cf.", "nr.", "aff.", "hybrid", "form", "n."]
	# changes items in removelist to 'sp.' if in blast definition
	defin = defstr.split()
	for x in removelist:
		if x in defin:
			defin[defin.index(x)] = "sp."
	# sometimes definitions look like 'PREDICTED: Genus Species'
	if defin[0] == "PREDICTED:":
		species = " ".join(defin[1:3])
	elif defin[0] != "PREDICTED:" and "." not in defin[0]:
		species = " ".join(defin[0:2])
	elif defin[0] != "PREDICTED:" and "." in defin[0]:
		species = defin[0].split(".")[0] + " " + defin[0].split(".")[1]
	definition = defstr
	# sometimes definitions have Lepidoptera sp. at the beginning and the species names are later - this parses these out
	# not generalized
	while "Lepidoptera" in species:
		try:
			defsplit = definition.split(">", 1)[-1]
			def2 = defsplit.split(' ', 1)[1]
			defin = def2.split()
			species = ' '.join(defin[0:2])
			definition = defsplit
		except:
			species = "Unknown species"
	# some more non generalized parsing - Burns and Robbins name things different from everything else
	if "sp." in species and "DHJ0" in defin[2] and "Burns0" not in defin[2] and "Robbins0" not in defin[2]:
		species = defin[0]+" "+defin[2].split("DHJ0")[0]
	elif "sp." in species and "DHJ0" in defin[2] and ("Burns0" in defin[2] or "Robbins0" in defin[2]):
		species = defin[0] + " " + defin[1]
	if "sp." in species and "ECO0" in defin[2]:
		species = defin[0]+ " " + defin[2].split("ECO0")[0]
	elif "sp." not in species and "ECO0" in defin[1]:
		species = defin[0] + " " + defin[1][:-5]
	if "BioLep" in defin[0]:
		species = defin[0][:-8] + " " + defin[1]

	genus = species.split()[0].replace("'", "")
	epithet = species.split()[1].replace("'", "")
	species = species.replace("'", "")
	return(species, genus, epithet)


