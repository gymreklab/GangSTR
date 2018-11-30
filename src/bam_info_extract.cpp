/*
Copyright (C) 2017 Melissa Gymrek <mgymrek@ucsd.edu>
and Nima Mousavi (mousavi@ucsd.edu)

This file is part of GangSTR.

GangSTR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GangSTR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GangSTR.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "src/bam_info_extract.h"
#include <iostream>
using namespace std;
BamInfoExtract::BamInfoExtract(const Options* options_,
			       BamCramMultiReader* bamreader_, 
			       RegionReader* region_reader_){
	options = options_;
	bamreader = bamreader_;
	region_reader = region_reader_;	
}

bool BamInfoExtract::GetReadLen(int32_t* read_len){
	*read_len = -1;
	int32_t flank_size = 20000;
	int32_t req_streak = 10;
	bool found_read_len = false, has_reads = false;
	// Header has info about chromosome names
	const BamHeader* bam_header = bamreader->bam_header();
	int32_t chrom_ref_id, num_reads;
	BamAlignment alignment;

	int32_t curr_len = 0, curr_streak = 0;

	while(region_reader->GetNextRegion(&locus) and !found_read_len){
		chrom_ref_id = bam_header->ref_id(locus.chrom);
		// if (chrom_ref_id == -1){
		// 	chrom_ref_id = bam_header->ref_id(locus.chrom.substr(3));
		// }

		// collecting reads mapped around locus
		bamreader->SetRegion(locus.chrom, 
			locus.start - flank_size > 0 ? locus.start - flank_size : 0, 
			locus.start + flank_size);
		num_reads = 0;
		curr_streak = 0;
		// Go through each alignment in the region until you have enough reads
		while (bamreader->GetNextAlignment(alignment) and curr_streak < req_streak) {
			has_reads = true;
			if(alignment.QueryBases().size() == curr_len){
				curr_streak++;
			}
			else{
				curr_len = alignment.QueryBases().size();
				curr_streak = 0;
			}
		}

		if (curr_streak >= req_streak){
			*read_len = curr_len;
			found_read_len = true;
		}
	}
	return found_read_len;
}

// TODO REWRITE THIS TO DEAL WITH PER-SAMPLE
// TODO Split to separate function to learn coverage, mean/sdev, and pdf/cdf
bool BamInfoExtract::GetInsertSizeDistribution(std::map<std::string, double>* sample_to_meandist,
					       std::map<std::string, double>* sample_to_sdev,
					       std::map<std::string, double>* sample_to_coverage,
					       std::map<std::string, std::vector<double> >* sample_to_pdf,
					       std::map<std::string, std::vector<double> >* sample_to_cdf,
					       const std::set<std::string> samples,
					       const std::map<std::string, std::string> rg_ids_to_sample) {
  double mean, std_dev, coverage;
  int32_t dist_size = options->dist_distribution_size;
  std::vector<double> dist_pdf(dist_size);
  std::vector<double> dist_cdf(dist_size);
  for (int i = 0; i < dist_size; i++){
    dist_pdf[i] = 0;
    dist_cdf[i] = 0;
  }
  // TODO change 200000 flank size to something appropriate
  int32_t flank_size = 400000;
  int32_t exclusion_margin = 1000;
  int32_t min_req_reads = 50000;
  bool found_ins_distribution = false, found_coverage = false;
  double mean_b, mean_a, std_b, std_a; // mean and std dev, before and after locus
  int* valid_temp_len_arr;
  std::vector<int32_t> temp_len_vec, valid_temp_len_vec;
  
  int reads_before = 0, reads_after = 0;
  // Header has info about chromosome names
  const BamHeader* bam_header = bamreader->bam_header();
  while(region_reader->GetNextRegion(&locus) and !found_ins_distribution and !found_coverage){
    const int32_t chrom_ref_id = bam_header->ref_id(locus.chrom);
    
    int32_t median, size = 0, sum = 0, valid_size = 0, sum_std = 0;
    BamAlignment alignment;
    
    // collecting reads mapped before locus
    bamreader->SetRegion(locus.chrom, 
			 locus.start - flank_size > 0 ? locus.start - flank_size : 0, 
			 locus.start - exclusion_margin > 0 ? locus.start - exclusion_margin : 0);
    
    // Go through each alignment in the region
    reads_before = 0;
    while (bamreader->GetNextAlignment(alignment)) {
      reads_before++;
      // Set template length
      temp_len_vec.push_back(abs(alignment.TemplateLength()));
      size++;
    }
    // collecting reads mapped after locus
    bamreader->SetRegion(locus.chrom, 
			 locus.start + exclusion_margin, 
			 locus.start + flank_size);
    // Go through each alignment in the region
    reads_after = 0;
    while (bamreader->GetNextAlignment(alignment)) {
      reads_after++;
      // Set template length
      temp_len_vec.push_back(abs(alignment.TemplateLength()));
      size++;
    }
    std::vector<int32_t> dist_count(dist_size);
    for (int i = 0; i < dist_size; i++){
      dist_count[i] = 0;
    }
    // if there's enough reads, compute and return TODO set threshold
    if (temp_len_vec.size() > min_req_reads) {
      if (reads_before > reads_after){
	coverage = float(reads_before * alignment.QueryBases().size()) / float(flank_size - exclusion_margin - alignment.QueryBases().size());
	found_coverage = true;
      }
      else {
	coverage = float(reads_after * alignment.QueryBases().size()) / float(flank_size - exclusion_margin - alignment.QueryBases().size());
	found_coverage = true;
      }
      
      sort(temp_len_vec.begin(), temp_len_vec.end());
      median = temp_len_vec.at(int32_t(size / 2));
      
      for (std::vector<int32_t>::iterator temp_it = temp_len_vec.begin();
	   temp_it != temp_len_vec.end();
	   ++temp_it) {
	// Todo change 3
	if(*temp_it < 4 * median and *temp_it > 0){
	  valid_temp_len_vec.push_back(*temp_it);
	  valid_size++;  
	  // Updating 
	  if (*temp_it < dist_size and *temp_it > 0){
	    dist_count[*temp_it]++;
	  }
	}
      }
      double cumulative = 0.0;
      dist_pdf[0] = double(dist_count[0] + dist_count[1] 
			   + dist_count[2]) / 5.0 / double(valid_size);
      cumulative += dist_pdf[0];
      dist_cdf[0] = cumulative;
      dist_pdf[1] = double(dist_count[0] + dist_count[1] 
			   + dist_count[2] + dist_count[3]) / 5.0 / double(valid_size);
      cumulative += dist_pdf[1];
      dist_cdf[1] = cumulative;
      ofstream ins_file;
      ins_file.open((options->outprefix + ".insdata.tab").c_str());
      for (int i = 2; i < dist_size - 2; i++){
	dist_pdf[i] = double(dist_count[i - 2] +
			     dist_count[i - 1] + 
			     dist_count[i] + 
			     dist_count[i + 1] + 
			     dist_count[i + 2]) / 5.0/ double(valid_size);
	cumulative += dist_pdf[i];
	dist_cdf[i] = cumulative;
	ins_file << i << "\t" << dist_pdf[i] << "\t" << dist_cdf[i] << endl;
      }
      
      dist_cdf[dist_size - 1] = 1.0;
      
      valid_temp_len_arr = &valid_temp_len_vec[0];
      mean = gsl_stats_int_mean(valid_temp_len_arr, 1, valid_size - 1);
      std_dev = gsl_stats_int_sd_m (valid_temp_len_arr,  1, valid_size, mean);
      found_ins_distribution = true;
    }
  }

  // TODO remove. Placeholder to set same value for all samples
  for (std::set<std::string>::const_iterator it=samples.begin();
       it != samples.end(); it++) {
    (*sample_to_meandist)[*it] = mean;
    (*sample_to_sdev)[*it] = std_dev;
    (*sample_to_coverage)[*it] = coverage;
    (*sample_to_pdf)[*it] = dist_pdf;
    (*sample_to_cdf)[*it] = dist_cdf;
  }
  return found_ins_distribution;
}

BamInfoExtract::~BamInfoExtract(){
}




