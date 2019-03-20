import sys

if len(sys.argv) < 3:
    print "Usage: python remove_messy.py in.bed out.bed supp1.bed,...,suppN.bed"
    sys.exit()
in_bed = sys.argv[1]
out_bed = sys.argv[2]
supp_bed = sys.argv[3].split(',')

disease_lines = {}

def is_in_between(line, next_line, target):
    # check matching chrom
    if line[0] != next_line[0] or line[0] != target[0] or next_line[0] != target[0]:
        return False
    # If target already exists, return False
    if line[1] == target[1] or next_line[1] == target[1]:
        return False
    if int(line[1]) < int(target[1]) and int(next_line[1]) > int(target[1]):
        return True
    return False

for bed in supp_bed:
    with open(bed, 'r') as disease_file:
        for rec in disease_file:
            recs = rec.strip().split('\t')
            line = [recs[0], recs[1], recs[2], \
                    recs[3], recs[4], \
                    ((int(recs[2]) - int(recs[1]) + 1) / int(recs[3])) * recs[4], str(0)]
            if recs[0] not in disease_lines:
                disease_lines[recs[0]] = []
            disease_lines[recs[0]].append(line)

        #print disease_lines

        
with open(in_bed, 'r') as refin:
    with open(out_bed, 'w') as refout:
        line = refin.readline().strip().split()
        while line != []:
            refout.write('\t'.join(line) + '\n')
            next_line = refin.readline().strip().split()
            if line[0] in disease_lines:
                for disease_line in disease_lines[line[0]]:
                    if next_line != [] and is_in_between(line, next_line, disease_line):
                        refout.write('\t'.join(disease_line) + '\n')        
            line = next_line
            
