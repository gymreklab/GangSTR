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

#include "src/options.h"
#include "src/region_reader.h"
#include "src/locus.h"
#include "src/bam_io.h"
#include "gsl/gsl_statistics_int.h"

#ifndef BAM_INFO_H_
#define BAM_INFO_H_

class BamInfoExtract{
public:
	BamInfoExtract(Options* options_,
						BamCramMultiReader* bamreader_, 
						RegionReader* region_reader_);
	~BamInfoExtract();
	bool GetReadLen(int32_t* read_len);
	bool GetInsertSizeDistribution(double* mean, double* std_dev, double *coverage);
private:
	Options* options;
	RegionReader* region_reader;
	Locus locus;
	BamCramMultiReader* bamreader;
};


#endif

