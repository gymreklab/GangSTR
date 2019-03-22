import argparse
import os,sys

def warning(s):
    print '>> WARNING: ' + s
def error(s):
    print '>> ERROR: ' + s
    sys.exit()

parser = argparse.ArgumentParser(description='Create STR-info file for GangSTR.')

parser.add_argument('--bed', type=str, required=True, help='GangSTR reference bed file.')
parser.add_argument('--out', type=str, required=True, help='Path to STR-info output')
parser.add_argument('--readlen', type=int, required=True, help='Read length')
parser.add_argument('--overwrite', type=str, help='A file in bed format containing thresholds that need to be overwritten')

args = parser.parse_args()
bed_file=args.bed
out_file=args.out
ov_bed=None
if not os.path.exists(bed_file):
    error('--bed file ' + bed_file + ' does not exist.')
out_dir='/'.join(out_file.split('/')[:-1])
if not os.path.exists(out_dir):
    error('--out dir ' + out_dir + ' does not exist.')
print out_file.split('/')[:-1]
asd
if args.overwrite is not None:
    ov_bed = args.overwrite
    if not os.path.exists(ov_bed):
        error('--overwrite bed file ' + ov_bed + ' does not exist.')
if args.readlen is not None:
    readlen = args.readlen
    if readlen <= 0:
        error('--readlen must be greater than 0.')

# Extract coordinates that need to be overwritten 
ov_dict={}
if ov_bed is not None:
    with open(ov_bed, 'r') as ov:
        i = 0
        for line in ov:
            i = i + 1
            recs = line.strip().split('\t')
            if len(recs) != 0 and len(recs) < 4:
                error ('Line ' + str(i) + ' of --overwrite bed does not have all the required columns,\
                or is not tab delimited: chrom start end threshold')
            chrom = recs[0]
            if chrom not in ov_dict:
                ov_dict[chrom] = {}
            pos_end = recs[1] + '_' + recs[2]
            if pos_end in ov_dict[chrom]:
                warning('Multiple lines in --overwrite bed for region ' + chrom + ':' + recs[1] + '-' + recs[2])
            ov_dict[chrom][pos_end] = recs[3]

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
            if chrom in ov_dict:
                if pos_end in ov_dict[chrom]:
                    thresh = ov_dict[chrom][pos_end]
                else:
                    thresh = int(float(readlen) / float(recs[3])) 
            else:
                thresh = int(float(readlen) / float(recs[3]))
            out_line = out_line + '\t' + thresh + '\n'
            out.write(out_line)
            





