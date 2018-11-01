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

#ifndef SRC_REALIGNMENT_H__
#define SRC_REALIGNMENT_H__

#include <string>
#include <vector>

#include <stdint.h>
#include "src/ssw_cpp.h"

// Set NW params
// NOTE: the other score triple (12,-12,-16) causes issues in realignment test (ambigous cases)
// Diag score too big ^^
// DO NOT USE! use SSW instead.
const static int32_t MATCH_SCORE = 3;
const static int32_t MISMATCH_SCORE = -1;
const static int32_t GAP_SCORE = -1;


// SSW Parameters
const static int32_t SSW_MATCH_SCORE = 2;
const static int32_t SSW_MISMATCH_SCORE = 5;
const static int32_t SSW_GAP_OPEN = 4;
const static int32_t SSW_GAP_EXTEND = 2;

//const static int32_t SSW_MATCH_SCORE = 2;
//const static int32_t SSW_MISMATCH_SCORE = 5;
//const static int32_t SSW_GAP_OPEN = 8;
//const static int32_t SSW_GAP_EXTEND = 2;


// amount of slip we allow between alignment position and STR start and end
static int32_t MARGIN = 5;		// This value is reset in expansion_aware_realign
// Threshold to discard alignment as non-overlapping
const static double MATCH_PERC_THRESHOLD = 0.88;

enum SingleReadType {
  SR_PREFLANK = 0,
  SR_POSTFLANK = 1,
  SR_ENCLOSING = 2,
  SR_IRR = 3,
  SR_UM_POT_IRR = 4, // Unmapped potential IRR
  SR_MAPPED_AFTER = 5,	// Not used TODO delete
  SR_UNKNOWN = 6,
  SR_NOT_FOUND = 7,		// Not used TODO delete
};

enum FlankMatchState{
  FM_NOMATCH = 0,
  FM_PARTIAL = 1,
  FM_COMPLETE = 2
};

enum sw_move{
	SW_END = 0,
	SW_DIAG = 1,
	SW_UP = 2,
	SW_LEFT = 3
};


bool trace_back(const std::vector<std::vector<int32_t> > score_matrix, 
            const int32_t& start_pos, 
            const int32_t& start_pos_temp,
            const std::string& seq1, 
            const std::string& seq2);

bool next_move(std::vector<std::vector<int32_t> > score_matrix, 
            const int32_t& x, 
            const int32_t& y, 
            sw_move* move);

bool expansion_aware_realign(const std::string& seq,
			     const std::string& qual,
			     const std::string& pre_flank,
			     const std::string& post_flank,
			     const std::string& motif,
			     const int32_t& min_match,
			     int32_t* nCopy, 
			     int32_t* start_pos, 
			     int32_t* end_pos, 
			     int32_t* score,
			     FlankMatchState* fm_start,
			     FlankMatchState* fm_end);		      

bool smith_waterman(const std::string& seq1,
		    const std::string& seq2,
		    const std::string& qual,
		    int32_t* pos, int32_t* end, int32_t* score);

bool striped_smith_waterman(const std::string& ref,
        const std::string& seq,
        const std::string& qual,
        int32_t* pos, int32_t* pos_temp, int32_t* score, int32_t* mismatches);
//static void ssw_PrintAlignment(const StripedSmithWaterman::Alignment& alignment);
bool create_score_matrix(const int32_t& rows, const int32_t& cols,
			 const std::string& seq1,
			 const std::string& seq2,
			 const std::string& qual,
			 std::vector<std::vector<int32_t> >* score_matrix,
			 int32_t* start_pos, int32_t* start_pos_temp, int32_t* current_score);

bool calc_score(const int32_t& i, const int32_t& j,
		const std::string& seq1, const std::string& seq2,
		const std::string& qual,
		std::vector<std::vector<int32_t> >* score_matrix);

bool classify_realigned_read(const std::string& seq,
			     const std::string& motif,
			     const int32_t& start_pos,
			     const int32_t& end_pos,
			     const int32_t& nCopy,
			     const int32_t& score,
			     const int32_t& prefix_length,
			     const int32_t& min_match,
			     const bool& isMapped,
			     const std::string& pre_flank,
			     const std::string& post_flank,
			     const FlankMatchState& fm_start,
			     const FlankMatchState& fm_end,
			     SingleReadType* single_read_class);

#endif  // SRC_REALIGNMENT_H__
