import sys

infile = sys.argv[1]
merge_thresh = int(sys.argv[2])
expansion = int(sys.argv[3])

def str_max (s1, s2):
    return str(max(int(s1), int(s2)))
def str_min (s1, s2):
    return str(min(int(s1), int(s2)))
def expand(loc):
    return [loc[0],\
            str(int(loc[1]) - expansion),\
            str(int(loc[2]) + expansion)]

def can_merge(loc1, loc2):
    if loc1[0] != loc2[0]:
        return False
    elif abs(int(loc1[1]) - int(loc2[1])) < merge_thresh or\
         abs(int(loc1[2]) - int(loc2[2])) < merge_thresh:
        return True
    return False
def merge(loc1, loc2):
    loc = [loc1[0],\
           str_min(loc1[1], loc2[1]),\
           str_max(loc1[2], loc2[2])]
    return loc

with open (infile, 'r') as file:
    loci = []
    for line in file:
        rec = line.strip().split('\t')
        merged = False
        for idx,loc in enumerate(loci):
            if can_merge(loci[idx], rec):
                #print loci[idx], rec
                loci[idx] = merge(loci[idx], rec)
                #print loci[idx]
                merged = True
                break
        if not merged:
            loci.append(rec)
            #print merged
            #print "Append", rec

out_string = ""
eh_string = ""
for idx, loc in enumerate(loci):
    if '_' in loc[0]:
        print 'Skipping ', loc
        continue
    eloc = expand(loc)
    out_string = out_string + eloc[0] + ':' +\
        eloc[1] + '-' +\
        eloc[2] + ','
print
print ">> Off target regions: "
print out_string[:-1]

