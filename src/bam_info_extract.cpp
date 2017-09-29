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
BamInfoExtract::BamInfoExtract(Options* options_,
						BamCramMultiReader* bamreader_, 
						RegionReader* region_reader_){
	options = options_;
	bamreader = bamreader_;
	region_reader = region_reader_;	
}

bool BamInfoExtract::GetReadLen(int32_t* read_len){
	int32_t flank_size = 2000;
	int32_t req_streak = 10;
	bool found_read_len = false;
	// Header has info about chromosome names
	const BamHeader* bam_header = bamreader->bam_header();
	int32_t chrom_ref_id, num_reads;
	BamAlignment alignment;

	int32_t curr_len, curr_streak = 0;

	while(region_reader->GetNextRegion(&locus) and !found_read_len){
		chrom_ref_id = bam_header->ref_id(locus.chrom);
		// collecting reads mapped around locus
		bamreader->SetRegion(locus.chrom, 
			locus.start - flank_size > 0 ? locus.start - flank_size : 0, 
			locus.start + flank_size);
		num_reads = 0;
		curr_streak = 0;
		// Go through each alignment in the region until you have enough reads
		while (bamreader->GetNextAlignment(alignment) and curr_streak < req_streak) {
			if(alignment.QueryBases().size() == curr_len){
				curr_streak++;
			}
			else{
				curr_len = alignment.QueryBases().size();
				curr_streak = 0;
			}
		}
		if (req_streak >= req_streak){
			*read_len = curr_len;
			found_read_len = true;
		}
	}
	if (found_read_len){
		return true;
	}
	else{
		*read_len = -1;
		return false;
	}
}

BamInfoExtract::~BamInfoExtract(){
}




