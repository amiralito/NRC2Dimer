#!/usr/bin/env python3

import sys
import math
from Bio import SeqIO
from collections import Counter


query_group_fasta = sys.argv[1]
subject_group_fasta = sys.argv[2]
min_freq_threshold_in_query = float(sys.argv[3])
min_ratio_threshold_in_query = float(sys.argv[4])

def check_thresholds(output_lines):
	for line in output_lines:
		freq = line.split(",")[1].split(":")[1]
		ratio = line.split(",")[2].split(":")[1]
		if freq != 'NA' and ratio != 'NA':
			if float(freq) > min_freq_threshold_in_query and float(ratio) > min_ratio_threshold_in_query:
				return True
	return False

query_group_records= SeqIO.parse(query_group_fasta, "fasta")
subject_group_records = SeqIO.parse(subject_group_fasta, "fasta")

query_group_seqs = [record.seq for record in query_group_records]
subject_group_seqs = [record.seq for record in subject_group_records]

total_n_seqs_in_query = len(query_group_seqs)
total_n_seqs_in_subject = len(subject_group_seqs)

total_alignment_length = len(query_group_seqs[0])

for n in range(total_alignment_length):
	query_group_aa_list_at_pos_n = [seq[n] for seq in query_group_seqs if seq[n] != "-"]
	subject_group_aa_list_at_pos_n = [seq[n] for seq in subject_group_seqs if seq[n] != "-"]
	query_group_aa_count_at_pos_n = Counter(query_group_aa_list_at_pos_n)
	subject_group_aa_count_at_pos_n = Counter(subject_group_aa_list_at_pos_n)
	output_lines = []
	if len(query_group_aa_count_at_pos_n) == 0:
		line = "-,Freq:NA,R_aa:NA"
		output_lines.append(line)
	else:
		for aa, query_cnt in sorted(query_group_aa_count_at_pos_n.items(), key=lambda x:x[1], reverse=True):
			query_aa_freq = query_cnt/total_n_seqs_in_query
			if aa not in subject_group_aa_count_at_pos_n.keys():
				line = "" + aa + ",Freq:" + str(round(query_aa_freq, 3)) + ',R_aa:inf'
			else:
				subject_aa_freq = subject_group_aa_count_at_pos_n[aa]/total_n_seqs_in_subject
				r_aa = round(query_aa_freq/subject_aa_freq, 3)
				line = "" + aa + ",Freq:" + str(round(query_aa_freq, 3)) + ',R_aa:' + str(r_aa)
			output_lines.append(line)
	
	if check_thresholds(output_lines):
		print("*", "Pos:" + str(n+1), ";".join(output_lines), sep="\t")
	else:
		print(".", "Pos:" + str(n+1), ";".join(output_lines), sep="\t")


