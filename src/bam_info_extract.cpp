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
BamInfoExtract::BamInfoExtract(Options* options_){
	options = options_;
	region_reader = new RegionReader(options->regionsfile);
	region_reader->GetNextRegion(&locus);
	int merge_type = BamCramMultiReader::ORDER_ALNS_BY_FILE;
  	bamreader = new BamCramMultiReader(options->bamfiles, options->reffa, merge_type);
	cerr<<locus.start;
}

int32_t BamInfoExtract::GetReadLen(){
	// TODO change 200000 flank size to something appropriate
	int32_t flank_size = 200000;
	int32_t exclusion_margin = 1000;

	// Header has info about chromosome names
	const BamHeader* bam_header = bamreader->bam_header();
	const int32_t chrom_ref_id = bam_header->ref_id(locus.chrom);
	BamAlignment alignment;

	collecting reads mapped before locus
	bamreader->SetRegion(locus.chrom, 
	  locus.start - flank_size > 0 ? locus.start - flank_size : 0, 
	  locus.start - exclusion_margin > 0 ? locus.start - exclusion_margin : 0);
	int32_t j = 0;
	Go through each alignment in the region
	while (bamreader->GetNextAlignment(alignment) and j < 100) {
		
		j++;
	}
}

BamInfoExtract::~BamInfoExtract(){
	delete region_reader;
	delete options;
	delete bamreader;
}




