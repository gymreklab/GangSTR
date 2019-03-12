import sys

if len(sys.argv) < 3:
    print "Usage: python minimal_trim.py in.bed out.bed"
    sys.exit()
in_bed = sys.argv[1]
out_bed = sys.argv[2]


# Check if the first and last K bps match, trim half motifs at the end.
def count_motif(repeat_str, motif):
    m = len(motif)
    c = 0
    s = 0
    while s < len(repeat_str):
        if (repeat_str[s:s + m] == motif):
            c = c + 1
            s = s + m
        else:
            s = s + 1
    return c


# Find coumpund motifs: ATAT = (AT)*2
def is_compound(motif):
    l = len(motif)
    threshold = 0.8
    for i in range (1, int(l / 2) + 1):
        sub = motif[0:i]
        # print count_motif(motif, sub) * i, l * threshold
        if count_motif(motif, sub) * i > l * threshold :
            return True
    return False


def minimal_trim(rep, motif):
    mm = motif * 2
    ll = len(motif) * 2
    start_match = False
    end_match = False
    max_trim_len = min(len(motif) * 3, int(len(rep) / 2))
    for start_offset in range(max_trim_len + 1):
        if (start_offset + ll >= len(rep)):
            return -1,-1
        if rep[start_offset: start_offset + ll] == mm:
            start_match = True
            break
    if start_match == False:
        return -1, -1
    for end_offset in list(reversed(range(len(rep) - max_trim_len, len(rep) + 1))):
        if end_offset - ll < 0:
            return -1, -1
        if rep[end_offset - ll: end_offset] == mm:
            end_match = True
            break
    if end_match == False:
        return -1, -1
    return start_offset, end_offset

thresholds = {1:10, 2:5, 3:4, 4:3, 5:3, 6:3}

with open(in_bed, 'r') as ref:
    with open(out_bed, 'w') as trim: # Save in temp to avoid overwriting   
        for line in ref:
            cols = line.strip().split('\t')

            ref_start = int(cols[1])
            ref_end = int(cols[2])
            motif = cols[4]
            m = len(motif)
            repeat_str = cols[5]

            if is_compound(motif):
                continue
            
            if len(motif) == 1:
                continue

            st,en = minimal_trim(repeat_str, motif)
            
            if st == -1 or en == -1:
                continue

            new_rep_str = repeat_str[st:en]
            cols[5] = new_rep_str
            cols[1] = str(ref_start + st)
            cols[2] = str(ref_start + en -1)
            ref_copy = (ref_end - ref_start + 1) / m


            if len(motif) in thresholds:
                thresh = thresholds[len(motif)]
            else:
                thresh = 3

            if (ref_copy >= thresh):
                out_line = "\t".join(cols) + "\t" + str(len(repeat_str) - len(new_rep_str)) + "\n"
                trim.write(out_line)

