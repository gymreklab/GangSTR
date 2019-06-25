import sys
import subprocess
help_dir = 'helpers/'
if len(sys.argv) < 4:
    print "Usage: python run.py motif ref_genome.fa temp_dir num_reads"
    sys.exit()

motif = sys.argv[1]
ref_genome = sys.argv[2]
temp_dir = sys.argv[3]
read_count = int(sys.argv[4])

read_len = 150
insert_size_mean = 500
insert_size_sdev = 100

helper_merge_thresh = 50
helper_expansion = 20

read_simulator = 'wgsim'
exp = 'ref_tmp'

# Create temp fasta
fa_path = temp_dir + exp + '.fa'
fa_len = 4000
with open(fa_path, 'w') as fa:
    fa.write('>' + motif + '\n')
    fa.write(motif * (fa_len / len(motif)))

# Simulate reads
fq1_path = temp_dir + exp + '.read1.fq'
fq2_path = temp_dir + exp + '.read2.fq'
base_error = 0.005

mutat_rate = 0.0
indel_frac = 0.0
indel_xtnd = 0.0


cmd = subprocess.Popen([read_simulator,
                        '-e', str(base_error),
                        '-d', str(insert_size_mean),
                        '-s', str(insert_size_sdev),
                        '-N', str(read_count),
                        '-1', str(read_len),
                        '-2', str(read_len),
                        '-r', str(mutat_rate),
                        '-R', str(indel_frac),
                        '-X', str(indel_xtnd),
                        '-S', str(0),
                        fa_path,
                        fq1_path,
                        fq2_path])
cmd.wait()


# Delete fasta
cmd = subprocess.Popen(['rm', fa_path])
cmd.wait()

# Align reads
sam_path = temp_dir + exp + '.sam'
sam_handle = open(sam_path, 'w')
cmd = subprocess.Popen(['bwa', 'mem',
                        '-M',
                        '-t', '4',
                        '-R', '@RG\\tID:foo\\tSM:bar',
                        ref_genome,
                        fq1_path,
                        fq2_path], stdout = sam_handle)
cmd.wait()
sam_handle.close()


# Delete fastq
cmd = subprocess.Popen(['rm', fq1_path])
cmd.wait()
cmd = subprocess.Popen(['rm', fq2_path])
cmd.wait()


# Process sam file
raw_regions = temp_dir + exp + '.txt'
cmd = subprocess.Popen(['bash', help_dir + 'process_sam.sh',\
                        sam_path,\
                        raw_regions])
cmd.wait()

# Delete sam
cmd = subprocess.Popen(['rm', sam_path])
cmd.wait()

# Merge regions
cmd = subprocess.Popen(['python', help_dir + 'merge_off.py', \
                        raw_regions, str(helper_merge_thresh), \
                        str(helper_expansion)])
cmd.wait()
# Delete raw regions
cmd = subprocess.Popen(['rm', raw_regions])
