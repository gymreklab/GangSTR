import subprocess, re
import numpy as np

locus = "HTT"
dataset_dir = "/storage/resources/datasets/repeat-expansions/bams/"
bed_dir = "../tests/" + locus + ".bed"
ref_genome = "/storage/resources/dbase/human/hg19/Homo_sapiens_assembly19.fasta"
bam_file_list = "bamlists/HTT_full.txt"

error_mode  = "rms"					# rms, mae: mean absolute error

sum_error_s = 0.0
sum_error_l = 0.0

samp_count = 0.0
with open (bam_file_list, 'r') as bam_list:
	for record in bam_list:
		bam_file = record.rstrip()

		real_genot = [float(x) for x in bam_file.split("/")[-1].split(".")[0].split("_")[1:3]]

		cmd = subprocess.Popen(["../src/GangSTR", \
						"--bam", bam_file, \
						"--ref", ref_genome, \
						"--regions", bed_dir, \
						"--out", "test"], stdout=subprocess.PIPE)
		cmd.wait()
		for line in cmd.stdout:
			result = line.rstrip()
			estm_genot = [float(x) for x in list(re.findall('(\d+)\, (\d+)',result))[0]]

		print '>>  ', (min(real_genot), max(real_genot)), '->\t', (min(estm_genot), max(estm_genot))
		samp_error_s = min(estm_genot) - min(real_genot)
		samp_error_l = max(estm_genot) - max(real_genot)
		samp_count = samp_count + 1.0

		if error_mode == "rms":
			sum_error_s = sum_error_s + samp_error_s * samp_error_s
			sum_error_l = sum_error_l + samp_error_l * samp_error_l
		elif error_mode == "mae":
			sum_error_s = sum_error_s + np.abs(samp_error_s)
			sum_error_l = sum_error_l + np.abs(samp_error_l)


		# if samp_count == 3:
		# 	break;

if error_mode == "rms":
	error_s = np.sqrt(sum_error_s / samp_count)
	error_l = np.sqrt(sum_error_l / samp_count)
	error = np.sqrt((sum_error_s + sum_error_l) / (2 * samp_count))
elif error_mode == "mae":
	error_s = sum_error_s / samp_count
	error_l = sum_error_l / samp_count
	error = sum_error_s + sum_error_l / (2 * samp_count)
print "Short allele " + error_mode + " error:\t", error_s
print "long allele " + error_mode + " error:\t", error_l
print "both alleles " + error_mode + " error:\t", error