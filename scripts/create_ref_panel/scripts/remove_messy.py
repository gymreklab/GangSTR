import ssw
import sys

if len(sys.argv) < 3:
    print "Usage: python remove_messy.py in.bed out.bed"
    sys.exit()
in_bed = sys.argv[1]
out_bed = sys.argv[2]


aligner = ssw.Aligner()
allowed_percent = 0.0
with open(in_bed, 'r') as refin:
    with open(out_bed, 'w') as refout:
        messy = 0
        clean = 0
        for line in refin:
            rec = line.strip().split("\t")
            ref = rec[5]
            query = rec[4] * int((int(rec[2]) - int(rec[1]) + 1) / int(rec[3]))            
            #alignment = aligner.align(ref, query)

            #if alignment.mismatch_count + alignment.deletion_count + alignment.insertion_count > allowed_percent / 100.0 * len(ref):
            if (ref != query):
                #print (alignment.alignment_report())
                messy = messy + 1
            else:
                clean = clean + 1
                refout.write('\t'.join(rec) + '\n')

print('Deleted ' + str(messy) + ' loci')
print(str(clean) + ' loci survived') 
            
    
