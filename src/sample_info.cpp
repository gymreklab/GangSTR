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

#include <iostream>

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
    PrintMessageDieOnError("Number of BAM files in --bam and samples in --bam-samps must match", M_ERROR, false);
  }
  for (size_t i=0; i<options.bamfiles.size(); i++) {
    PrintMessageDieOnError("Loading read group  " + read_groups[i] + " for file " + options.bamfiles[i], M_PROGRESS, options.quiet);
    rg_ids_to_sample[options.bamfiles[i]] = read_groups[i];
    rg_samples.insert(read_groups[i]);
  }
  return true;
}

bool SampleInfo::LoadReadGroups(const Options& options, const BamCramMultiReader& bamreader) {
  for (size_t i=0; i<options.bamfiles.size(); i++) {
    const std::vector<ReadGroup>& read_groups = bamreader.bam_header(i)->read_groups();
    if (read_groups.empty()) {
      PrintMessageDieOnError("\tNo read group specified in BAM file", M_ERROR, false);
    } 
    for (std::vector<ReadGroup>::const_iterator rg_iter = read_groups.begin(); rg_iter != read_groups.end(); rg_iter++) {
      if (!rg_iter->HasID()) {
	PrintMessageDieOnError("RG in BAM/CRAM header is lacking the ID tag", M_ERROR, false);
      }
      if (!rg_iter->HasSample()) {
	PrintMessageDieOnError("RG in BAM/CRAM header is lacking the SM tag",M_ERROR, false);
      }
      if (rg_ids_to_sample.find(rg_iter->GetID()) != rg_ids_to_sample.end()) {
	if (rg_ids_to_sample[rg_iter->GetID()].compare(rg_iter->GetSample()) != 0) {
	  PrintMessageDieOnError("Read group id " + rg_iter->GetID() + " maps to more than one sample", M_ERROR, false);
	}
      }
      PrintMessageDieOnError("Loading read group id " + rg_iter->GetID() + " for sample " + rg_iter->GetSample(), M_PROGRESS, options.quiet);
      rg_ids_to_sample[options.bamfiles[i]+":"+rg_iter->GetID()] = rg_iter->GetSample();
      rg_samples.insert(rg_iter->GetSample());
    } 
  }
  return true;
}

bool SampleInfo::ExtractBamInfo(const Options& options, BamCramMultiReader& bamreader,
				RegionReader& region_reader, const RefGenome& ref_genome) {
  BamInfoExtract bam_info(&options, &bamreader, &region_reader, &ref_genome);
  if (options.read_len == -1) {
    if (!bam_info.GetReadLen(&read_len)) {
      PrintMessageDieOnError("Error extracting read length", M_ERROR, false);
    }
  }
  else{
    read_len = options.read_len;
  }
  region_reader.Reset();
  
  // Set insert size distribution
  if (options.dist_mean.empty() or options.dist_sdev.empty()) {
    if (!bam_info.GetInsertSizeDistribution(&profile, rg_samples, rg_ids_to_sample, custom_read_groups)) {
      PrintMessageDieOnError("Error extracting insert size info", M_ERROR, false);
    }
  } else {
    if (options.dist_mean.size() != options.dist_sdev.size()) {
      PrintMessageDieOnError("Different number of dist means and dist sdevs input", M_ERROR, false);
    }
    if (options.dist_mean.size() == 1) {
      size_t i = 0;
      for (std::set<std::string>::iterator it = rg_samples.begin();
	   it != rg_samples.end(); it++) {
	profile[*it].dist_mean = options.dist_mean[i];
	profile[*it].dist_sdev = options.dist_sdev[i];
	// Precalculating pdf and cdf
	int32_t dist_size = options.dist_distribution_size;
	std::vector<double> dist_pdf(dist_size);
	std::vector<double> dist_cdf(dist_size);
	std::vector<double> dist_integral(dist_size);
	double integral = 0;
	for (int j = 0; j < dist_size; j++){
	  dist_pdf[j] = gsl_ran_gaussian_pdf(j - options.dist_mean[i], options.dist_sdev[i]);
	  dist_cdf[j] = gsl_cdf_gaussian_P(j - options.dist_mean[i], options.dist_sdev[i]);
	  integral += j * dist_pdf[j];
	  dist_integral[j] = integral;
	}
	dist_cdf[dist_size - 1] = 1.0;
	profile[*it].dist_pdf = dist_pdf;
	profile[*it].dist_cdf = dist_cdf;
	profile[*it].dist_integral = dist_integral;
      }
    } else {
      if (options.dist_mean.size() != options.bamfiles.size()) {
	PrintMessageDieOnError("Different number of dist means and BAM files input", M_ERROR, false);
      }
      if (!custom_read_groups) {
	PrintMessageDieOnError("Can only set per-BAM dists if using custom read groups", M_ERROR, false);
      }
      for (size_t i=0; i<options.bamfiles.size(); i++) {
	cerr << "set for: " << options.bamfiles[i] << endl;
	// Precalculating pdf and cdf
	int32_t dist_size = options.dist_distribution_size;
	std::vector<double> dist_pdf(dist_size);
	std::vector<double> dist_cdf(dist_size);
	std::vector<double> dist_integral(dist_size);
	double integral = 0;
	for (int j = 0; j < dist_size; j++){
	  dist_pdf[j] = gsl_ran_gaussian_pdf(j - options.dist_mean[i], options.dist_sdev[i]);
	  dist_cdf[j] = gsl_cdf_gaussian_P(j - options.dist_mean[i], options.dist_sdev[i]);
	  integral += j * dist_pdf[j];
	  dist_integral[j] = integral;
	}
	dist_cdf[dist_size - 1] = 1.0;
	profile[GetSampleFromID(options.bamfiles[i])].dist_mean = options.dist_mean[i];
	profile[GetSampleFromID(options.bamfiles[i])].dist_sdev = options.dist_sdev[i];
	profile[GetSampleFromID(options.bamfiles[i])].dist_pdf = dist_pdf;
	profile[GetSampleFromID(options.bamfiles[i])].dist_cdf = dist_cdf;
	profile[GetSampleFromID(options.bamfiles[i])].dist_integral = dist_integral;
      }
    }
  }

  // Set coverage
  if (options.coverage.empty()) {
    if (!bam_info.GetCoverage(&profile, rg_samples, rg_ids_to_sample, custom_read_groups)) {
      PrintMessageDieOnError("Error extracting coverage info", M_ERROR, false);
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
	PrintMessageDieOnError("Different number of coverages and BAM files input", M_ERROR, false);
      }
      if (!custom_read_groups) {
	PrintMessageDieOnError("Can only set per-BAM coverages if using custom read groups", M_ERROR, false);
      }
      for (size_t i=0; i<options.bamfiles.size(); i++) {
	profile[options.bamfiles[i]].coverage = options.coverage[i];
      }
    }
  }
  return true;
}

void SampleInfo::PrintSampleInfo(const std::string& logfilename) {
  stringstream ss;
  int numgc;
  ss << "Sample stats:\n";
  for (std::set<std::string>::iterator it = rg_samples.begin();
       it != rg_samples.end(); it++) {
    ss << *it << "\n"
       << "\tCoverage=" << profile[*it].coverage << "\n"
       << "\tInsMean=" << profile[*it].dist_mean << "\n"
       << "\tInsSdev=" << profile[*it].dist_sdev << "\n"
       <<"\tReadLen="<< read_len << "\n";
    numgc = profile[*it].gc_coverage.size();
    if (numgc > 0) {
      ss << "\tGC Coverage" << "\n";
      for (size_t i=0; i<profile[*it].gc_coverage.size();i++) {
	ss << "\t\t Bin " << i << " " << profile[*it].gc_coverage[i] << "\n";
      }
    }
  }
  //  PrintMessageDieOnError(ss.str(), M_PROGRESS,);
  // now print to log file
  ofstream sample_log_file;
  sample_log_file.open(logfilename.c_str());
  sample_log_file << "Sample\tMeanCoverage\tInsMean\tInsSdev\tReadLen";
  for (int i=0; i<numgc; i++) {
    sample_log_file << "\tCoverage-GCBin-"<< i;
  }
  sample_log_file << "\n";
  for (std::set<std::string>::iterator it=rg_samples.begin();
       it != rg_samples.end(); it++) {
    sample_log_file << *it << "\t"
	     << profile[*it].coverage << "\t"
	     << profile[*it].dist_mean << "\t"
	     << profile[*it].dist_sdev << "\t"
	     << read_len;
    for (int i=0; i<numgc; i++) {
      sample_log_file << "\t" << profile[*it].gc_coverage[i];
    }
    sample_log_file << "\n";
  }
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

const double SampleInfo::GetGCCoverage(std::string sample, int32_t gcbin) {
  if (gcbin == -1 || gcbin >= profile[sample].gc_coverage.size()) {
    //    PrintMessageDieOnError("GC Bin not set or out of range. Returning mean coverage", M_WARNING);
    return profile[sample].coverage;
  }
  return profile[sample].gc_coverage[gcbin];
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

bool SampleInfo::GetSampleProfile(const std::string& sample, SampleProfile* sp) {
  if (profile.find(sample) != profile.end()) {
    *sp = profile[sample];
    return true;
  }
  return false;
}

SampleInfo::~SampleInfo() {
}
