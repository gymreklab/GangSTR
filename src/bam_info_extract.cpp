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
#include "src/gc_region_reader.h"
#include <iostream>
using namespace std;
BamInfoExtract::BamInfoExtract(const Options* options_,
			       BamCramMultiReader* bamreader_, 
			       RegionReader* region_reader_,
			       const RefGenome* ref_genome_){
	options = options_;
	bamreader = bamreader_;
	region_reader = region_reader_;
	ref_genome = ref_genome_;
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

bool BamInfoExtract::GetCoverageGC(std::map<std::string, SampleProfile>* profile,
				   const std::set<std::string> samples,
				   std::map<std::string, std::string> rg_ids_to_sample,
				   bool custom_read_groups) {
  int regions_per_bin = 1000;
  float lb, ub;
  std::vector<int32_t> total_bases;
  std::map<std::string, std::vector<int32_t> > sample_gc_bases;
  GCRegionReader gc_reader(*ref_genome, options->gc_bin_size, options->gc_region_len,
			   options->max_gc_regions);
  std::string read_group, rgid, sample, fname;
  BamAlignment alignment;
  bool found_sample;
  for (int i=0; i<(1.0/options->gc_bin_size); i++) {
    // Init each sample for this bin
    total_bases.push_back(0);
    for (std::set<std::string>::const_iterator sampleit=samples.begin(); sampleit != samples.end(); sampleit++) {
      sample_gc_bases[*sampleit].push_back(0);
    }
    // Get regions
    lb = options->gc_bin_size*i;
    ub = options->gc_bin_size*(i+1);
    std::vector<Locus> gc_bin_loci;
    if (!gc_reader.GetGCBinLoci(&gc_bin_loci, lb, ub,
				regions_per_bin)) {
      stringstream ss;
      ss << "Could not find GC regions " << lb << "-" << ub;
      PrintMessageDieOnError(ss.str(), M_WARNING, options->quiet);
    }
    // For each locus get num total bases and coverage per sample
    for (std::vector<Locus>::iterator it = gc_bin_loci.begin(); it != gc_bin_loci.end(); it++) {
      const Locus locus = *it;
      bamreader->SetRegion(locus.chrom, locus.start, locus.end);
      //      std::cerr << "Locus " << locus.chrom << ":" << locus.start << " " << lb <<"-" << ub << std::endl;
      while (bamreader->GetNextAlignment(alignment)) {
	// Is this read worth looking at?
	if (alignment.IsSupplementary() || alignment.IsSecondary() || 
	    !alignment.IsProperPair()) {
	  continue;
	}
	// What sample did this read come from?
	fname = alignment.file_;
	if (custom_read_groups) {
	  rgid = fname;
	  found_sample = true;
	} else {
	  if (!alignment.GetStringTag("RG", read_group)) {
	    PrintMessageDieOnError("Could not find read group for " + alignment.Name(), M_WARNING, options->quiet);
	    found_sample = false;
	  }
	  rgid = fname + ":" + read_group;
	  found_sample = true;
	}
	if (!found_sample) continue;
	if (rg_ids_to_sample.find(rgid) == rg_ids_to_sample.end()) continue;
	sample = rg_ids_to_sample[rgid];
	// Get template length and assign to that sample
	sample_gc_bases[sample][i] += alignment.QueryBases().size();
      }
      total_bases[i] += options->gc_region_len;
    }
  }

  for (std::set<std::string>::const_iterator it = samples.begin();
       it != samples.end(); it++) {
    std::vector<double> sample_gc_covs;
    for (int i = 0; i<total_bases.size(); i++) {
      if (total_bases[i]>0) {
	sample_gc_covs.push_back(float(sample_gc_bases[*it][i])/float(total_bases[i]));
      } else {
	sample_gc_covs.push_back(-1);
      }
    }
    (*profile)[*it].gc_coverage = sample_gc_covs;
  }
}

bool BamInfoExtract::GetCoverage(std::map<std::string, SampleProfile>* profile,
				 const std::set<std::string> samples,
				 std::map<std::string, std::string> rg_ids_to_sample,
				 bool custom_read_groups) {
  if (options->model_gc_cov) {
    GetCoverageGC(profile, samples, rg_ids_to_sample, custom_read_groups);
  }
  // How many regions etc. to use
  int num_regions_to_use = 100;
  int num_regions_so_far = 0;
  int region_offset = 10000; // Look this far away from STR
  int region_length = 5000; // Use this length of region to look at
  int total_region_length = 0;
  // Keep track of coverage for each sample
  std::map<std::string, int> sample_to_bases;
  for (std::set<std::string>::const_iterator it = samples.begin();
       it != samples.end(); it++) {
    sample_to_bases[*it] = 0;
  }
  // Set up
  std::string read_group, rgid, sample, fname;
  bool found_sample;
  Locus locus;
  region_reader->Reset();
  BamAlignment alignment;
  while ((num_regions_so_far < num_regions_to_use) &&
	 region_reader->GetNextRegion(&locus)) {
    bamreader->SetRegion(locus.chrom, locus.start+region_offset, locus.start+region_offset+region_length);
    while (bamreader->GetNextAlignment(alignment)) {
      // Is this read worth looking at?
      if (alignment.IsSupplementary() || alignment.IsSecondary() || 
	  !alignment.IsProperPair()) {
	continue;
      }
      // What sample did this read come from?
      fname = alignment.file_;
      if (custom_read_groups) {
	rgid = fname;
	found_sample = true;
      } else {
	if (!alignment.GetStringTag("RG", read_group)) {
	  PrintMessageDieOnError("Could not find read group for " + alignment.Name(), M_WARNING, options->quiet);
	  found_sample = false;
	}
	rgid = fname + ":" + read_group;
	found_sample = true;
      }
      if (!found_sample) continue;
      if (rg_ids_to_sample.find(rgid) == rg_ids_to_sample.end()) continue;
      sample = rg_ids_to_sample[rgid];
      // Get template length and assign to that sample
      sample_to_bases[sample] += alignment.QueryBases().size();
    }
    num_regions_so_far++;
    total_region_length += region_length;
  }
  for (std::set<std::string>::const_iterator it = samples.begin();
       it != samples.end(); it++) {
    (*profile)[*it].coverage = float(sample_to_bases[*it])/float(total_region_length);
  }

  return true;
}

bool BamInfoExtract::GetInsertSizeDistribution(std::map<std::string, SampleProfile>* profile,
					       const std::set<std::string> samples,
					       std::map<std::string, std::string> rg_ids_to_sample,
					       bool custom_read_groups) {
  // Keep track of template lengths for each sample

  std::map<std::string, std::vector<int32_t> > sample_to_tlens;
  for (std::set<std::string>::const_iterator it = samples.begin();
       it != samples.end(); it++) {
    std::vector<int32_t> vec;
    vec.clear();
    sample_to_tlens[*it] = vec;
  }

  // How many regions etc. to use
  int num_regions_to_use = 100;
  int num_regions_so_far = 0;
  int region_offset = 1000; // Look this far away from STR
  int region_length = 5000; // Use this length of region to look at
  // Requirements to continue with each sample
  size_t min_reads_per_sample = options->min_reads_per_sample / 2 * (options->ploidy==-1?2:options->ploidy);
  // Set up
  std::string read_group, rgid, sample, fname;
  bool found_sample;
  Locus locus;
  region_reader->Reset();
  BamAlignment alignment;
  while ((num_regions_so_far < num_regions_to_use) &&
	 region_reader->GetNextRegion(&locus)) {
    bamreader->SetRegion(locus.chrom, locus.start+region_offset, locus.start+region_offset+region_length);
    while (bamreader->GetNextAlignment(alignment)) {
      // Is this read worth looking at?
      if (alignment.IsSupplementary() || alignment.IsSecondary() || 
	  !alignment.IsProperPair()) {
	continue;
      }
      // What sample did this read come from?
      fname = alignment.file_;
      if (custom_read_groups) {
	rgid = fname;
	found_sample = true;
      } else {
	if (!alignment.GetStringTag("RG", read_group)) {
	  PrintMessageDieOnError("Could not find read group for " + alignment.Name(), M_WARNING, options->quiet);
	  found_sample = false;
	}
	rgid = fname + ":" + read_group;
	found_sample = true;
      }
      if (!found_sample) continue;
      if (rg_ids_to_sample.find(rgid) == rg_ids_to_sample.end()) continue;
      sample = rg_ids_to_sample[rgid];
      // Get template length and assign to that sample
      sample_to_tlens[sample].push_back(abs(alignment.TemplateLength()));
    }
    num_regions_so_far++;
  }
  
  // Summarize distributions
  ofstream ins_file;
  ins_file.open((options->outprefix + ".insdata.tab").c_str());  
  for (std::set<std::string>::const_iterator it=samples.begin();
       it != samples.end(); it++) {
    std::vector<int32_t> tlen_vec = sample_to_tlens[*it];
    // Set up to get values
    int32_t dist_size = options->dist_distribution_size;
    std::vector<int32_t> dist_count(dist_size);
    std::vector<double> dist_pdf(dist_size);
    std::vector<double> dist_cdf(dist_size);
    std::vector<double> dist_integral(dist_size);
    for (int i = 0; i < dist_size; i++){
      dist_count[i] = 0;
      dist_pdf[i] = 0;
      dist_cdf[i] = 0;
      dist_integral[i] = 0;
    }

    if (tlen_vec.size() < min_reads_per_sample) {
      std::stringstream ss;
      ss << "Not enough reads for " << *it << " " << tlen_vec.size()<<". Please set insert size distribution manually.";
      PrintMessageDieOnError(ss.str(), M_ERROR, false);
      return false;
    }
    // Filter out extremes
    sort(tlen_vec.begin(), tlen_vec.end());
    int32_t median = tlen_vec.at(int32_t(tlen_vec.size()/2));
    std::vector<int32_t> tlen_vec_filt;
    for (size_t i=0; i<tlen_vec.size(); i++) {
      if (tlen_vec[i]>0 && tlen_vec[i]<4*median) {
	tlen_vec_filt.push_back(tlen_vec[i]);
	if (tlen_vec[i]<dist_size) {
	  dist_count[tlen_vec[i]]++;
	}
      }
    }
    if (tlen_vec_filt.size() < min_reads_per_sample) {
      std::stringstream ss;
      ss << "Not enough reads for " << *it << " " << tlen_vec_filt.size()<<". Please set insert size distribution manually.";
      PrintMessageDieOnError(ss.str(), M_ERROR, false);
      return false;
    }
    // Set pdf/cdf
    double cumulative = 0.0;
    double integral = 0.0;
    double smoother = 5.0; // Smooth together this many bins
    double total_reads = double(tlen_vec_filt.size());
    dist_pdf[0] = double(dist_count[0] + dist_count[1] 
			 + dist_count[2]) / smoother / total_reads;
    cumulative += dist_pdf[0];
    dist_cdf[0] = cumulative;
    dist_integral[0] = integral;
    dist_pdf[1] = double(dist_count[0] + dist_count[1] 
			 + dist_count[2] + dist_count[3]) / smoother / total_reads;
    cumulative += dist_pdf[1];
    integral += 1 * dist_pdf[1];
    dist_cdf[1] = cumulative;
    dist_integral[0] = integral;
    for (int i = 2; i < dist_size - 2; i++){
      dist_pdf[i] = double(dist_count[i - 2] +
			   dist_count[i - 1] + 
			   dist_count[i] + 
			   dist_count[i + 1] + 
			   dist_count[i + 2]) / smoother / total_reads;
      cumulative += dist_pdf[i];
      integral += i * dist_pdf[i];
      dist_cdf[i] = cumulative;
      dist_integral[i] = integral;
      ins_file << *it << " " << i << "\t" << dist_pdf[i] << "\t" << dist_cdf[i] << endl;
    }
      
    dist_cdf[dist_size - 1] = 1.0;
      
    // Set values in profile
    (*profile)[*it].dist_mean = gsl_stats_int_mean(&tlen_vec_filt[0], 1, tlen_vec_filt.size());
    (*profile)[*it].dist_sdev = gsl_stats_int_sd_m(&tlen_vec_filt[0], 1, tlen_vec_filt.size(),
						   (*profile)[*it].dist_mean);
    (*profile)[*it].dist_pdf = dist_pdf;
    (*profile)[*it].dist_cdf = dist_cdf;
    (*profile)[*it].dist_integral = dist_integral;
  }
  return true;
}

BamInfoExtract::~BamInfoExtract(){
}




