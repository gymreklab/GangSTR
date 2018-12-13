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

#ifndef SRC_SAMPLE_INFO_H__
#define SRC_SAMPLE_INFO_H__

#include <map>
#include <set>
#include <string>
#include <vector>

#include "src/bam_info_extract.h"
#include "src/options.h"
#include "src/region_reader.h"


/*
  Keep track of per-sample info on read groups,
  template length distributions, and coverage
 */
class SampleInfo {
  const static int dist_size = 2000;
 public:
  SampleInfo();
  virtual ~SampleInfo();

  /* Functions for setting sample info */
  bool SetCustomReadGroups(const Options& options);
  bool LoadReadGroups(const Options& options, const BamCramMultiReader& bamreader);
  bool ExtractBamInfo(const Options& options, BamCramMultiReader& bamreader,
		      RegionReader& region_reader, const RefGenome& ref_genome);

  /* Getters */
  const int32_t GetReadLength();
  const std::set<std::string> GetSamples();
  const double GetInsertMean(std::string sample);
  const double GetInsertSdev(std::string sample);
  const double GetCoverage(std::string sample);
  const bool GetIsCustomRG();
  std::vector<double> GetDistPDF(std::string sample);
  std::vector<double> GetDistCDF(std::string sample);
  std::vector<double> GetDistIntegral(std::string sample);
  const std::string GetSampleFromID(const std::string& rgid);

  /* Other utils */
  void PrintSampleInfo(const std::string& logfilename);
  double GetDistMax(const std::string& sample);

 private:
  bool custom_read_groups;
  std::set<std::string> rg_samples;
  std::map<std::string, std::string> rg_ids_to_sample;
  std::map<std::string, SampleProfile> profile;
  
  int32_t read_len;
};

#endif  // SRC_SAMPLE_INFO_H__
