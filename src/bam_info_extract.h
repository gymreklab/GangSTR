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

#include <set>

#include "src/options.h"
#include "src/region_reader.h"
#include "src/locus.h"
#include "src/bam_io.h"
#include "gsl/gsl_statistics_int.h"

#ifndef BAM_INFO_H_
#define BAM_INFO_H_

class BamInfoExtract{
public:
	BamInfoExtract(const Options* options_,
		       BamCramMultiReader* bamreader_, 
		       RegionReader* region_reader_);
	~BamInfoExtract();
	bool GetReadLen(int32_t* read_len);
	// TODO change to deal with per sample
	bool GetInsertSizeDistribution(std::map<std::string, double>* sample_to_meandist,
				       std::map<std::string, double>* sample_to_sdev,
				       std::map<std::string, double>* sample_to_coverage,
				       std::map<std::string, double*>* sample_to_pdf,
				       std::map<std::string, double*>* sample_to_cdf,
				       const std::set<std::string> samples,
				       const std::map<std::string, std::string> rg_ids_to_sample);
 private:
	const Options* options;
	BamCramMultiReader* bamreader;
	RegionReader* region_reader;
	Locus locus;
};


#endif

