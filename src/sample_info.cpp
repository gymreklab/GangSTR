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
  else{
    read_len = options.read_len;
  }
  region_reader.Reset();
  
  // Set insert size distribution
  if (options.dist_mean.empty() or options.dist_sdev.empty()) {
    if (!bam_info.GetInsertSizeDistribution(&profile, rg_samples, rg_ids_to_sample, custom_read_groups)) {
      PrintMessageDieOnError("Error extracting insert size info", M_ERROR);
    }
  } else {
    if (options.dist_mean.size() != options.dist_sdev.size()) {
      PrintMessageDieOnError("Different number of dist means and dist sdevs input", M_ERROR);
    }
    if (options.dist_mean.size() == 1) {
      size_t i = 0;
      for (std::set<std::string>::iterator it = rg_samples.begin();
	   it != rg_samples.end(); it++) {
	profile[*it].dist_mean = options.dist_mean[i];
	profile[*it].dist_sdev = options.dist_sdev[i];
      }
    } else {
      if (options.dist_mean.size() != options.bamfiles.size()) {
	PrintMessageDieOnError("Different number of dist means and BAM files input", M_ERROR);
      }
      if (!custom_read_groups) {
	PrintMessageDieOnError("Can only set per-BAM dists if using custom read groups", M_ERROR);
      }
      for (size_t i=0; i<options.bamfiles.size(); i++) {
	profile[options.bamfiles[i]].dist_mean = options.dist_mean[i];
	profile[options.bamfiles[i]].dist_sdev = options.dist_sdev[i];
      }
    }
  }

  // Set coverage
  if (options.coverage.empty()) {
    if (!bam_info.GetCoverage(&profile, rg_samples, rg_ids_to_sample, custom_read_groups)) {
      PrintMessageDieOnError("Error extracting coverage info", M_ERROR);
    }
  } else {
    if (options.coverage.size() == 1) {
      size_t i = 0;
      for (std::set<std::string>::iterator it = rg_samples.begin();
	   it != rg_samples.end(); it++) {
	profile[*it].coverage = options.coverage[i];
      }
    } else {
      if (options.coverage.size() != options.bamfiles.size()) {
	PrintMessageDieOnError("Different number of coverages and BAM files input", M_ERROR);
      }
      if (!custom_read_groups) {
	PrintMessageDieOnError("Can only set per-BAM coverages if using custom read groups", M_ERROR);
      }
      for (size_t i=0; i<options.bamfiles.size(); i++) {
	profile[options.bamfiles[i]].coverage = options.coverage[i];
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
       << "\tCoverage=" << profile[*it].coverage << "\n"
       << "\tInsMean=" << profile[*it].dist_mean << "\n"
       << "\tInsSdev=" << profile[*it].dist_sdev << "\n"
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
  return profile[sample].dist_mean;
}

const double SampleInfo::GetInsertSdev(std::string sample) {
  return profile[sample].dist_sdev;
}

const double SampleInfo::GetCoverage(std::string sample) {
  return profile[sample].coverage;
}

std::vector<double> SampleInfo::GetDistPDF(std::string sample) {
  return profile[sample].dist_pdf;
}

std::vector<double> SampleInfo::GetDistCDF(std::string sample) {
  return profile[sample].dist_cdf;
}

std::vector<double> SampleInfo::GetDistIntegral(std::string sample) {
  return profile[sample].dist_integral;
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
  dist_mean = profile[sample].dist_mean;
  dist_sdev = profile[sample].dist_sdev;
  return dist_mean+3*dist_sdev;
}

SampleInfo::~SampleInfo() {
}
