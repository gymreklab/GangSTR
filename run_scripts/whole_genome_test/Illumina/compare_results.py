import subprocess

bam_dir = "/storage/resources/datasets/IlluminaRepeatExpansions/"
repo_dir = "/storage/nmmsv/GangSTR"
ref_genome = "/storage/resources/dbase/human/hg19/hg19.fa"

genotype={}
with open('true_genotypes.txt', 'r') as genots_file:
	for line in genots_file:
		cols = line.strip().split('\t')
		genotype[cols[0]] = (cols[1], cols[4])

bam={}
for sample in genotype:
	ls = subprocess.Popen(['ls', bam_dir], stdout=subprocess.PIPE)
	bam_name = subprocess.Popen(['grep', sample], stdin=ls.stdout, stdout=subprocess.PIPE)
	for line in bam_name.stdout:
		str_line = line.strip()
		if str_line[-3:] == 'bai':
			bam[sample] = str_line[:-4]
			print bam[sample]
		# print sample, line.strip().split('.cip')[0]
	# print sample, genotype[sample]

for sample in bam:
	locus = genotype[sample][0]
	grnd_truth = genotype[sample][1]
	bam_file = bam_dir + bam[sample]
	print bam[sample].split('_')[4][:-4]	
	print grnd_truth
	cmd = subprocess.Popen([repo_dir + "/src/GangSTR", \
								"--bam", bam_file, \
								"--ref", ref_genome, \
								"--regions", 'loci/'+locus+'.bed', \
								"--out", "test", \
								"--frrweight", str(0.5), \
								"--enclweight", "1.0", \
								"--spanweight", "1.0", \
								"--flankweight", str(1.0),\
								"--ploidy", str(2),\
								"--numbstrap", str(100),\
								"--minmatch", str(4),\
								"--minscore", str(80)])
	cmd.wait()
	print "###"
	
