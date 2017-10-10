import subprocess, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# locus = "HTT"
locus = "CACNA1A"
# locus = "ATXN7"
repo_dir = "/storage/nmmsv/GangSTR"
run_dir = repo_dir + "/run_scripts/"
# result_dir = run_dir + "/results2/"
# dataset_dir = "/storage/resources/datasets/repeat-expansions/bams/"
bed_dir = repo_dir + "/tests/" + locus + ".bed"
ref_genome = "/storage/resources/dbase/human/hg19/Homo_sapiens_assembly19.fasta"
# bam_file_list = run_dir + "bamlists/HTT_full.txt"
bam_file_list = run_dir + "bamlists/CACNA1A.txt"
# bam_file_list = run_dir + "bamlists/ATXN7_grid.txt"
true_available = True

error_mode  = "rms"					# rms, mae: mean absolute error

cmd = subprocess.Popen(["make", "-j", "-C", repo_dir])
cmd.wait()


sum_error_s = 0.0
sum_error_l = 0.0
samp_count = 0.0

error_s_list = []
error_l_list = []
with open (bam_file_list, 'r') as bam_list:
	with open('genot_errors.txt', 'w') as results:
		for record in bam_list:
			bam_file = record.rstrip()

			if true_available == True:
				real_genot = [float(x) for x in bam_file.split("/")[-1].split(".")[0].split("_")[1:3]]

			cmd = subprocess.Popen([repo_dir + "/src/GangSTR", \
									"--bam", bam_file, \
									"--ref", ref_genome, \
									"--regions", bed_dir, \
									"--out", "test", \
									"--frrweight", str(0.5), \
									"--enclweight", "1.0", \
									"--spanweight", "1.0", \
									"--flankweight", str(1.0),\
									"--ploidy", str(2),\
									"--numbstrap", str(100),\
									"--minmatch", str(4),\
									"--minscore", str(80)], stdout=subprocess.PIPE)
			cmd.wait()
			for line in cmd.stdout:
				result = line.rstrip()
				if line[0] == ">":
					estm_genot = [float(x) for x in list(re.findall('(\d+)\, (\d+)',result))[0]]
				else:
					print result
			if true_available == False:				
				print '>>  ', (min(estm_genot), max(estm_genot))
			else:
				print '>>  ', (min(real_genot), max(real_genot)), '->\t', (min(estm_genot), max(estm_genot))
				samp_error_s = min(estm_genot) - min(real_genot)
				samp_error_l = max(estm_genot) - max(real_genot)
				error_s_list.append(samp_error_s)
				error_l_list.append(samp_error_l)
				samp_count = samp_count + 1.0

				if error_mode == "rms":
					sum_error_s = sum_error_s + samp_error_s * samp_error_s
					sum_error_l = sum_error_l + samp_error_l * samp_error_l
				elif error_mode == "mae":
					sum_error_s = sum_error_s + np.abs(samp_error_s)
					sum_error_l = sum_error_l + np.abs(samp_error_l)

				print samp_error_s, samp_error_l, samp_count
			
			results.write(bam_file.split("/")[-1] + "\t" +\
							str(int(min(estm_genot))) + "\t" +\
							str(int(max(estm_genot))) + "\n")
			# if samp_count == 15:
			# 	break;
print
print ">> Result for", locus

if true_available == True:
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

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.hist(error_s_list, normed=False, bins=30)
	ax.set_ylabel('Absolute Genotyping Error')
	ax.set_title('Shorter Allele')
	fig.savefig('genot_error_hist_s.pdf')

	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111)
	ax2.hist(error_l_list, normed=False, bins=20)
	ax2.set_ylabel('Absolute Genotyping Error')
	ax2.set_title('Longer Allele')
	fig2.savefig('genot_error_hist_l.pdf')