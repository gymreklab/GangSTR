import csv

with open('true_genotypes.txt', 'w') as output:
	with open('Supplemental_Table_7.csv', 'r') as csvfile:
		for record in csvfile:
			lines = record.split('\r')
			for line in lines:
				cols = line.split(',')
				if cols[0] != '':
					if cols[0][0] == 'N' or cols[0][0] == 'C':
						rec = '\t'.join([cols[0], cols[1], cols[2], cols[3], cols[7]]) + '\n'
						output.write(rec)
				

	

