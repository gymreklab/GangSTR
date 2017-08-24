import subprocess, re, os, errno
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# StackOverflow
def mkdir_p(path):
    try:
        os.makedirs(path)
        print 'Created path: ', path
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


locus = "HTT"
repo_dir = "/storage/nmmsv/GangSTR"
run_dir = repo_dir + "/run_scripts/"
result_dir = run_dir + "/results3_after_fixing_extract/"
dataset_dir = "/storage/resources/datasets/repeat-expansions/bams/"
bed_dir = repo_dir + "/tests/" + locus + ".bed"
ref_genome = "/storage/resources/dbase/human/hg19/Homo_sapiens_assembly19.fasta"
bam_file_list = run_dir + "bamlists/HTT_full.txt"
true_available = True

error_mode  = "rms"					# rms, mae: mean absolute error

cmd = subprocess.Popen(["make", "-j", "-C", repo_dir])
cmd.wait()


limit = 11.0
list1 = np.arange(0,0.3,0.05)
list0 = np.arange(0,0.3,0.05)
heat_array_short = np.zeros((len(list0), len(list1)))
heat_array_long = np.zeros((len(list0), len(list1)))

i = 0
for frr_w in list0:
	j = 0
	for flank_w in list1:
		sum_error_s = 0.0
		sum_error_l = 0.0
		samp_count = 0
		error_s_list = []
		error_l_list = []
		exp_dir = result_dir + "frr" + str(int(frr_w * 100)) + \
				"_flank" + str(int(flank_w * 100))
		mkdir_p(exp_dir)
		print
		print "## Running for frr: " + str(frr_w) + ", flank: " + str(flank_w)
		print 
		with open (bam_file_list, 'r') as bam_list:
			with open(exp_dir + '/genotypes.txt', 'w') as results:
				for record in bam_list:
					bam_file = record.rstrip()

					if true_available == True:
						real_genot = [float(x) for x in bam_file.split("/")[-1].split(".")[0].split("_")[1:3]]

					cmd = subprocess.Popen([repo_dir + "/src/GangSTR", \
									"--bam", bam_file, \
									"--ref", ref_genome, \
									"--regions", bed_dir, \
									"--out", "test", \
									"--frrweight", str(frr_w), \
									"--enclweight", "0.3", \
									"--spanweight", "1.0", \
									"--flankweight", str(flank_w)], stdout=subprocess.PIPE)
					cmd.wait()

					for line in cmd.stdout:
						result = line.rstrip()
						estm_genot = [float(x) for x in list(re.findall('(\d+)\, (\d+)',result))[0]]
					if true_available == False:				
						print '>>  ', (min(estm_genot), max(estm_genot))
					else:
						print '>>  ', (min(real_genot), max(real_genot)), '->\t', (min(estm_genot), max(estm_genot))
						samp_error_s = min(estm_genot) - min(real_genot)
						samp_error_l = max(estm_genot) - max(real_genot)
						error_s_list.append(samp_error_s)
						error_l_list.append(samp_error_l)
						
						

						samp_count = samp_count + 1

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
					# if samp_count == 3:
					# 	break;

		if true_available == True:
			if error_mode == "rms":
				error_s = np.sqrt(sum_error_s / samp_count)
				error_l = np.sqrt(sum_error_l / samp_count)
				error = np.sqrt((sum_error_s + sum_error_l) / (2 * samp_count))
			elif error_mode == "mae":
				error_s = sum_error_s / samp_count
				error_l = sum_error_l / samp_count
				error = sum_error_s + sum_error_l / (2 * samp_count)

			with open(exp_dir + "/error.txt", "w") as error_file:
				error_file.write("Short\t" + error_mode + "\t" + str(error_s) + "\n")
				error_file.write("Long\t" + error_mode + "\t" + str(error_l) + "\n")
				error_file.write("Both\t" + error_mode + "\t" + str(error))

			print "Short allele " + error_mode + " error:\t", error_s
			print "long allele " + error_mode + " error:\t", error_l
			print "both alleles " + error_mode + " error:\t", error

			if error_s < limit:
				heat_array_short[i][j] = error_s
			else:
				heat_array_short[i][j] = limit

			if error_l < limit:
				heat_array_long[i][j] = error_l
			else:
				heat_array_long[i][j] = limit

			fig = plt.figure()
			ax = fig.add_subplot(111)
			ax.hist(error_s_list, normed=False, bins=30)
			ax.set_ylabel('Absolute Genotyping Error')
			ax.set_title('Shorter Allele')
			fig.savefig(exp_dir + '/genot_error_hist_s.pdf')

			fig2 = plt.figure()
			ax2 = fig2.add_subplot(111)
			ax2.hist(error_l_list, normed=False, bins=30)
			ax2.set_ylabel('Absolute Genotyping Error')
			ax2.set_title('Longer Allele')
			fig2.savefig(exp_dir + '/genot_error_hist_l.pdf')
		j = j + 1
	i = i + 1

fig3 = plt.figure()
ax = fig3.add_subplot(111)
cax = ax.imshow(heat_array_short, cmap = 'hot_r', interpolation = 'nearest')
ax.set_ylabel('FRR Weight')
ax.set_xlabel('Flanking Weight')
ax.set_title('Shorter Allele | Error Heat Map')

plt.yticks(np.arange(len(list0)))
labels = [item.get_text() for item in ax.get_yticklabels()]
labels = [str(item) for item in list0]
ax.set_yticklabels(labels)

plt.xticks(np.arange(len(list1)))
labels = [item.get_text() for item in ax.get_xticklabels()]
labels = [str(item) for item in list1]
ax.set_xticklabels(labels)

fig3.colorbar(cax)
fig3.savefig(result_dir + "/short.pdf")

fig4 = plt.figure()
ax = fig4.add_subplot(111)
cax = ax.imshow(heat_array_long, cmap = 'hot_r', interpolation = 'nearest')
ax.set_ylabel('FRR Weight')
ax.set_xlabel('Flanking Weight')
ax.set_title('Longer Allele | Error Heat Map')

plt.yticks(np.arange(len(list0)))
labels = [item.get_text() for item in ax.get_yticklabels()]
labels = [str(item) for item in list0]
ax.set_yticklabels(labels)

plt.xticks(np.arange(len(list1)))
labels = [item.get_text() for item in ax.get_xticklabels()]
labels = [str(item) for item in list1]
ax.set_xticklabels(labels)


fig4.colorbar(cax)
fig4.savefig(result_dir + "/long.pdf")


