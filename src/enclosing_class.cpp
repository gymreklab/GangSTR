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

#include "src/mathops.h"
#include "src/enclosing_class.h"

#include <math.h>

using namespace std;

bool EnclosingClass::GetLogClassProb(const int32_t& allele,
				     double* log_class_prob) {
	int flank_len = 2000;	// arg_dict['flank_len']			// (F)
	int read_len =  100;		// arg_dict['read_len']				// (r)
	int motif_len = 3;		// arg_dict['motif']
	int str_len = allele * motif_len;					// (L)

	if (read_len <= str_len)
		*log_class_prob = 0;
	else
		*log_class_prob = 2.0 * double(read_len - str_len) / double(2 * flank_len + str_len - 2 * read_len);
 	return false; // TODO
}

bool EnclosingClass::GetLogReadProb(const int32_t& allele,
				    const int32_t& data,
				    double* log_allele_prob) {
	// FILLLLLL INNNNN Stutter model
	double u = 0.01;
	double d = 0.02;
	double p = 0.95;
	double delta = data - allele;
	if (delta == 0)
		*log_allele_prob = 1 - u - d;
	else if (delta > 0)
		*log_allele_prob = pow(u * p * (1.0 - p), (delta - 1.0));
	else
		*log_allele_prob = pow(d * p * (1.0 - p), (-delta - 1.0));
  	return false; // TODO
}

