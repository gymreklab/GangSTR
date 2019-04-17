import argparse
import os,sys

def warning(s):
    print '>> WARNING: ' + s
def error(s):
    print '>> ERROR: ' + s
    sys.exit()
nucToNumber={"A":0,"C":1,"G":2,"T":3}
def getCanonicalMS(repseq):
    """ Get canonical STR sequence """
    size = len(repseq)
    canonical = repseq
    for i in range(size):
        newseq = repseq[size-i:]+repseq[0:size-i]
        for j in range(size):
            if nucToNumber[newseq[j]] < nucToNumber[canonical[j]]:
                canonical = newseq
            elif nucToNumber[newseq[j]] > nucToNumber[canonical[j]]:
                break
    return canonical

parser = argparse.ArgumentParser(description='Create STR-info file for GangSTR.')

parser.add_argument('--bed', type=str, required=True, help='GangSTR reference bed file.')
parser.add_argument('--out', type=str, required=True, help='Path to STR-info output')
parser.add_argument('--readlen', type=int, required=True, help='Read length used to calculate threshold (Priority 1)')
parser.add_argument('--population', type=str, help='A file in bed format containing thresholdsf for low priority regions that need to be overwritten (Priority 2)')
parser.add_argument('--motif', nargs=2, action='append', help='Overwrite all canonical instances of the motif to a value. format: --motif motif1 thresh1 --motif motif2 thresh2 (Priority e)')
parser.add_argument('--disease', type=str, help='A file in bed format containing thresholds for high priority regions that need to be overwritten (Priority 4)')
parser.add_argument('--population-pad', type=int, default=0, help='Incresing the thresholds in population bed by this amount to avoind false positives.')





args = parser.parse_args()

bed_file=args.bed
out_file=args.out
ov_bed=None
pop_bed=None
if not os.path.exists(bed_file):
    error('--bed file ' + bed_file + ' does not exist.')
out_dir='/'.join(out_file.split('/')[:-1])
if not os.path.exists(out_dir):
    error('--out dir ' + out_dir + ' does not exist.')


if args.disease is not None:
    ov_bed = args.disease
    if not os.path.exists(ov_bed):
        error('--disease bed file ' + ov_bed + ' does not exist.')
if args.population is not None:
    pop_bed = args.population
    if not os.path.exists(pop_bed):
        error('--population bed file ' + pop_bed + ' does not exist.')
    pop_pad = args.population_pad
if args.readlen is not None:
    readlen = args.readlen
    if readlen <= 0:
        error('--readlen must be greater than 0.')
motif_dict = {}
if args.motif is not None:
    raw_motif_dict = dict(args.motif)
    for key in raw_motif_dict:
        try:
            num_thresh = int(raw_motif_dict[key])
        except:
            error ('Threshold for motif ' + key + ' is not a number: ' +raw_motif_dict[key])
        if num_thresh <= 0:
            error ('Threshold for motif ' + key + ' set to non-positive value: ' +raw_motif_dict[key])
        else:
            motif_dict[getCanonicalMS(key)] = num_thresh


# Extract coordinates that need to be overwritten 
ov_dict={}
if ov_bed is not None:
    with open(ov_bed, 'r') as ov:
        i = 0
        for line in ov:
            i = i + 1
            recs = line.strip().split('\t')
            # Skip header
            if recs[0]=='chrom':
                continue
            if len(recs) != 0 and len(recs) < 4:
                error ('Line ' + str(i) + ' of --disease bed does not have all the required columns,\
                or is not tab delimited: chrom start end threshold')
            chrom = recs[0]
            if chrom not in ov_dict:
                ov_dict[chrom] = {}
            pos_end = recs[1] + '_' + recs[2]
            if pos_end in ov_dict[chrom]:
                warning('Multiple lines in --disease bed for region ' + chrom + ':' + recs[1] + '-' + recs[2])
            ov_dict[chrom][pos_end] = recs[3]
pop_dict={}
if pop_bed is not None:
    with open(pop_bed, 'r') as pop:
        i = 0
        for line in pop:
            i = i + 1
            recs = line.strip().split('\t')
            # Skip header
            if recs[0]=='chrom':
                continue
            if len(recs) != 0 and len(recs) < 4:
                error ('Line ' + str(i) + ' of --population bed does not have all the required columns,\
                or is not tab delimited: chrom start end threshold')
            chrom = recs[0]
            if chrom not in pop_dict:
                pop_dict[chrom] = {}
            pos_end = recs[1] + '_' + recs[2]
            if pos_end in pop_dict[chrom]:
                warning('Multiple lines in --population bed for region ' + chrom + ':' + recs[1] + '-' + recs[2])
            pop_dict[chrom][pos_end] = str(int(recs[3]) + pop_pad)
header = '\t'.join(['chrom','pos','end','thresh'])
# Open bed file and create str-info for each line
with open(out_file, 'w') as out:
    out.write(header + '\n')
    with open(bed_file, 'r') as bed:
        j = 0
        for line in bed:
            j = j + 1
            recs = line.strip().split('\t')
            if len(recs) != 0 and len(recs) < 5:
                error ('Line ' + str(j) + ' of --bed does not have all the required columns,\
                or is not tab delimited: chrom start end motif_len motif')
            out_line = '\t'.join(recs[0:-2])
            chrom = recs[0]
            pos_end = recs[1] + '_' + recs[2]
            motif = recs[4]
            canonical_motif = getCanonicalMS(motif)
            # Priority 1 (lowest): Set thresh based on read length
            thresh = int(float(readlen) / float(recs[3])) 
            # Priority 2: Set thresh based on population bed
            if chrom in pop_dict:
                if pos_end in pop_dict[chrom]:
                    thresh = pop_dict[chrom][pos_end]
            # Priority 3: Set thresh based on motif
            if canonical_motif in motif_dict:
                thresh = motif_dict[canonical_motif]
            # Priority 4 (highest): Set thresh based on disease bed
            if chrom in ov_dict:
                if pos_end in ov_dict[chrom]:
                    thresh = ov_dict[chrom][pos_end]

            out_line = out_line + '\t' + str(thresh) + '\n'
            out.write(out_line)
            





