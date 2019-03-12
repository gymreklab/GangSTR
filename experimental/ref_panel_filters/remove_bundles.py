import sys

if len(sys.argv) < 4:
    print "Usage: python remove_bundles.py in.bed out.bed threshold"
    sys.exit()
in_bed = sys.argv[1]
out_bed = sys.argv[2]
thresh = int(sys.argv[3]) #50

discard_line = False
check_motif = False # False for new, True for loose


def is_close(line1, line2, check_motif):
    if check_motif and line1[4] != line2[4]:
        return False
    if line1[0] != line2[0]:
        return False
    s1 = int(line1[1])
    s2 = int(line2[1])
    e1 = int(line1[2])
    e2 = int(line2[2])
    if abs(s1 - s2) < thresh or\
       abs(s1 - e2) < thresh or\
       abs(s2 - e1) < thresh or\
       abs(e2 - e1) < thresh:
        return True
    
    return False


with open(in_bed, 'r') as refin:
    with open(out_bed, 'w') as refout:
        line = refin.readline().strip().split()
        while line != []:
            next_line = refin.readline().strip().split()
            if next_line != [] and is_close(line, next_line, check_motif):
                discard_line = True
                while (next_line != [] and is_close(line, next_line, check_motif)):
                    line = next_line
                    next_line = refin.readline().strip().split()
            #print line
            if not discard_line:
                refout.write('\t'.join(line) + '\n')
            discard_line = False
            line = next_line
            
    
#print "NO NEED TO RUN DEDUP WHEN USING THIS SCRIPT"
