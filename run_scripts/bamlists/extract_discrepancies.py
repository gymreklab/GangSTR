with open ('HTT_difficult.txt', 'w') as f_out:
	with open('eh_results_compare.txt', 'r') as f_in:
		for record in f_in:
			col = record.split('\t')
			if col[1] == col[2] and col[3] != col[1]:	# Expansion hunter correct, GangSTR wrong
				f_out.write('/storage/resources/datasets/repeat-expansions/bams/' + col[0] + '.rmdup.bam\n')
