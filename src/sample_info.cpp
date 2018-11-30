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

#include "src/sample_info.h"
#include "src/stringops.h"

using namespace std;

SampleInfo::SampleInfo() {
  custom_read_groups = false;
  rg_samples.clear();
  rg_ids_to_sample.clear();
  
}
/*
SampleInfo::SampleInfo(SampleInfo samp) {
  custom_read_groups = samp.custom_read_groups;
  rg_samples = samp.rg_samples;
  rg_ids_to_sample = samp.rg_ids_to_sample;
  sample_to_meandist = samp.sample_to_meandist;
  sample_to_sdev = samp.sample_to_sdev;
  sample_to_coverage = samp.sample_to_coverage;
  sample_to_pdf = new
}
*/
bool SampleInfo::SetCustomReadGroups(const Options& options) {
  custom_read_groups = true;
  std::vector<std::string> read_groups;
  split_by_delim(options.rg_sample_string, ',', read_groups);

  if (options.bamfiles.size() != read_groups.size()) {
    PrintMessageDieOnError("Number of BAM files in --bam and samples in --bam-samps must match", M_ERROR);
  }
  for (size_t i=0; i<options.bamfiles.size(); i++) {
    PrintMessageDieOnError("Loading read group  " + read_groups[i] + " for file " + options.bamfiles[i], M_PROGRESS);
    rg_ids_to_sample[options.bamfiles[i]] = read_groups[i];
    rg_samples.insert(read_groups[i]);
  }
  return true;
}

bool SampleInfo::LoadReadGroups(const Options& options, const BamCramMultiReader& bamreader) {
  for (size_t i=0; i<options.bamfiles.size(); i++) {
    const std::vector<ReadGroup>& read_groups = bamreader.bam_header(i)->read_groups();
    if (read_groups.empty()) {
      PrintMessageDieOnError("\tNo read group specified in BAM file", M_ERROR);
    } 
    for (std::vector<ReadGroup>::const_iterator rg_iter = read_groups.begin(); rg_iter != read_groups.end(); rg_iter++) {
      if (!rg_iter->HasID()) {
	PrintMessageDieOnError("RG in BAM/CRAM header is lacking the ID tag", M_ERROR);
      }
      if (!rg_iter->HasSample()) {
	PrintMessageDieOnError("RG in BAM/CRAM header is lacking the SM tag",M_ERROR);
      }
      if (rg_ids_to_sample.find(rg_iter->GetID()) != rg_ids_to_sample.end()) {
	if (rg_ids_to_sample[rg_iter->GetID()].compare(rg_iter->GetSample()) != 0) {
	  PrintMessageDieOnError("Read group id " + rg_iter->GetID() + " maps to more than one sample", M_ERROR);
	}
      }
      PrintMessageDieOnError("Loading read group id " + rg_iter->GetID() + " for sample " + rg_iter->GetSample(), M_PROGRESS);
      rg_ids_to_sample[options.bamfiles[i]+":"+rg_iter->GetID()] = rg_iter->GetSample();
      rg_samples.insert(rg_iter->GetSample());
    } 
  }
  return true;
}

bool SampleInfo::ExtractBamInfo(const Options& options, BamCramMultiReader& bamreader,
				RegionReader& region_reader) {
  BamInfoExtract bam_info(&options, &bamreader, &region_reader);
  if (options.read_len == -1) {
    if (!bam_info.GetReadLen(&read_len)) {
      PrintMessageDieOnError("Error extracting read length", M_ERROR);
    }
  }
  region_reader.Reset();
  // Insert size distribution + coverage inferred per sample
  if(options.dist_mean.empty() or options.dist_sdev.empty() or options.coverage.empty()){
    if(!bam_info.GetInsertSizeDistribution(&sample_to_meandist, &sample_to_sdev, &sample_to_coverage,
					   &sample_to_pdf, &sample_to_cdf, rg_samples, rg_ids_to_sample)) {
      PrintMessageDieOnError("Error extracting insert size and coverage info", M_ERROR);
    }
  }

  // Deal with setting custom ins/coverages per BAM file
  if (!options.dist_mean.empty() & !options.dist_sdev.empty()) {
    if (options.dist_mean.size() != options.dist_sdev.size()) {
      PrintMessageDieOnError("Different number of dist means and dist sdevs input", M_ERROR);
    }
    if (options.dist_mean.size() == 1) {
      size_t i = 0;
      for (std::set<std::string>::iterator it = rg_samples.begin();
	   it != rg_samples.end(); it++) {
	sample_to_meandist[*it] = options.dist_mean[i];
	sample_to_sdev[*it] = options.dist_sdev[i];
      }
    } else {
      if (options.dist_mean.size() != options.bamfiles.size()) {
	PrintMessageDieOnError("Different number of dist means and BAM files input", M_ERROR);
      }
      if (!custom_read_groups) {
	PrintMessageDieOnError("Can only set per-BAM dists if using custom read groups", M_ERROR);
      }
      for (size_t i=0; i<options.bamfiles.size(); i++) {
	sample_to_meandist[options.bamfiles[i]] = options.dist_mean[i];
	sample_to_sdev[options.bamfiles[i]] = options.dist_sdev[i];
      }
    }
  }
  if (options.use_cov & !options.coverage.empty()) {
    if (options.coverage.size() == 1) {
      size_t i = 0;
      for (std::set<std::string>::iterator it = rg_samples.begin();
	   it != rg_samples.end(); it++) {
	sample_to_coverage[*it] = options.coverage[i];
      }
    } else {
      if (options.coverage.size() != options.bamfiles.size()) {
	PrintMessageDieOnError("Different number of coverages and BAM files input", M_ERROR);
      }
      if (!custom_read_groups) {
	PrintMessageDieOnError("Can only set per-BAM coverages if using custom read groups", M_ERROR);
      }
      for (size_t i=0; i<options.bamfiles.size(); i++) {
	sample_to_coverage[options.bamfiles[i]] = options.coverage[i];
      }
    }
  }
  return true;
}

void SampleInfo::PrintSampleInfo() {
  stringstream ss;
  ss << "Sample stats:\n";
  for (std::set<std::string>::iterator it = rg_samples.begin();
       it != rg_samples.end(); it++) {
    ss << *it << "\n"
       << "\tCoverage=" << sample_to_coverage[*it] << "\n"
       << "\tInsMean=" << sample_to_meandist[*it] << "\n"
       << "\tInsSdev=" << sample_to_sdev[*it] << "\n"
       <<"\tReadLen="<< read_len << "\n";
  }
  PrintMessageDieOnError(ss.str(), M_PROGRESS);
}

const int32_t SampleInfo::GetReadLength() {
  return read_len;
}

const std::set<std::string> SampleInfo::GetSamples() {
  return rg_samples;
}

const double SampleInfo::GetInsertMean(std::string sample) {
  for (map<std::string,double>::iterator it = sample_to_meandist.begin(); 
       it != sample_to_meandist.end(); ++it){
    //cerr << it->first << endl;
  } 
  return sample_to_meandist[sample];
}

const double SampleInfo::GetInsertSdev(std::string sample) {
  return sample_to_sdev[sample];
}

const double SampleInfo::GetCoverage(std::string sample) {
  return sample_to_coverage[sample];
}

std::vector<double> SampleInfo::GetDistPDF(std::string sample) {
  return sample_to_pdf[sample];
}

std::vector<double> SampleInfo::GetDistCDF(std::string sample) {
  return sample_to_cdf[sample];
}

const bool SampleInfo::GetIsCustomRG() {
  return custom_read_groups;
}

const std::string SampleInfo::GetSampleFromID(const std::string& rgid) {
  for (map<std::string,std::string>::iterator it = rg_ids_to_sample.begin(); 
       it != rg_ids_to_sample.end(); ++it){
    //cerr << it->first << " " << it->second<< endl;
  } 
  return rg_ids_to_sample[rgid];
}

double SampleInfo::GetDistMax(const std::string& sample) {
  double dist_mean, dist_sdev;
  dist_mean = sample_to_meandist[sample];
  dist_sdev = sample_to_sdev[sample];
  return dist_mean+3*dist_sdev;
}

SampleInfo::~SampleInfo() {
}
