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

#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "src/str_info.h"
#include "src/common.h"
#include "src/stringops.h"

using namespace std;

STRInfo::STRInfo(const Options& options) {
  default_info.exp_thresh = 40;
  default_info.stutter_up = options.stutter_up;
  default_info.stutter_down = options.stutter_down;
  default_info.stutter_p = options.stutter_p;

  // Load params from file
  if (!options.str_info_file.empty()) {
     std::ifstream freader(options.str_info_file.c_str());
     if (!freader.is_open()) {
       PrintMessageDieOnError("Could not open STR info file" , M_ERROR);
     }
     // Get header and col nums
     std::string line;
     std::vector<std::string> items;
     int chrom_col, start_col, end_col, thresh_col,
       stutter_up_col, stutter_down_col, stutter_p_col;
     if (!std::getline(freader, line)) {
       PrintMessageDieOnError("No header in STR info file" , M_ERROR);
     }
     split_by_delim(line, '\t', items);
     if (items.size() <= 3) {
       PrintMessageDieOnError("No relevant STR info found. Check format." , M_ERROR);
     }
     chrom_col = 0;
     start_col = 1;
     thresh_col = -1;
     stutter_up_col = -1;
     stutter_down_col = -1;
     stutter_p_col = -1;
     for (size_t i=3; i<items.size(); i++) {
       const std::string colname = items[i];
       if (colname == "thresh") {
	 thresh_col = i;
       } else if (colname == "stutter_up") {
	 stutter_up_col = i;
       } else if (colname == "stutter_down") {
	 stutter_down_col = i;
       } else if (colname == "stutter_p") {
	 stutter_p_col = i;
       } else {
	 PrintMessageDieOnError("Unknown STR info column detected... " + colname, M_WARNING);
       }
     }
     // Get info for each STR
     std::string chrom;
     int32_t start, thresh;
     double stutter_up, stutter_down, stutter_p;
     while (std::getline(freader, line)) {
       items.clear();
       split_by_delim(line, '\t', items);
       chrom = items[chrom_col];
       start = atoi(items[start_col].c_str());
       STRLocusInfo sli;
       sli.exp_thresh = default_info.exp_thresh;
       sli.stutter_up = default_info.stutter_up;
       sli.stutter_down = default_info.stutter_down;
       sli.stutter_p = default_info.stutter_p;
       if (thresh_col != -1) {
	 sli.exp_thresh = atof(items[thresh_col].c_str());
       }
       if (stutter_up_col != -1) {
	 sli.stutter_up = atof(items[stutter_up_col].c_str());
       }
       if (stutter_down_col != -1) {
	 sli.stutter_up = atof(items[stutter_down_col].c_str());
       }
       if (stutter_p_col != -1) {
	 sli.stutter_p = atof(items[stutter_p_col].c_str());
       }
       std::pair<std::string, int32_t> locus(chrom, start);
       str_info[locus] = sli;
     }
  }
}

const int32_t STRInfo::GetExpansionThreshold(const std::string& chrom, const int32_t& start) {
  std::pair<std::string, int32_t> locus(chrom, start);
  if (str_info.find(locus) == str_info.end()) {
    return default_info.exp_thresh;
  } else {
    return str_info[locus].exp_thresh;
  }
}

STRInfo::~STRInfo() {}
