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
#include <iostream>
#include <algorithm>
#include <map>
using namespace std;

bool EnclosingClass::GetLogClassProb(const int32_t& allele,
				     const int32_t& read_len, const int32_t& motif_len,
				     double* log_class_prob) {
	int str_len = allele * motif_len;					// (L)
	double class_prob;
	if (read_len <= str_len)
		class_prob = 0;
	else
		class_prob = double(read_len - str_len) / double(2 * flank_len + str_len - 2 * read_len);

	if (class_prob > 0){
		*log_class_prob = log(class_prob);
		return true;
	}
	else if (class_prob == 0){
		*log_class_prob = NEG_INF;
		return true;
	}
	else
		return false;
}

bool EnclosingClass::GetLogReadProb(const int32_t& allele,
				    const int32_t& data,
				    const int32_t& read_len,
				    const int32_t& motif_len,
				    const int32_t& ref_count,
				    double* log_allele_prob) {
	double delta = data - allele;
	double allele_prob;
	if (delta == 0)
		allele_prob = 1 - stutter_up - stutter_down;
	else if (delta > 0)
		allele_prob = stutter_up * stutter_p * pow(1.0 - stutter_p, delta - 1.0);
	else
		allele_prob = stutter_down * stutter_p * pow(1.0 - stutter_p, -delta - 1.0);
	if (allele_prob > 0){
		*log_allele_prob = log(allele_prob);
		return true;
	}
	else if (allele_prob == 0){
		*log_allele_prob = NEG_INF;
		return true;
	}
	else
		return false;
}

bool EnclosingClass::ExtractEnclosingAlleles(std::vector<int> *alleles){
    std::map<int32_t, int32_t> allele_repeats;

	for (std::vector<int32_t>::iterator data_it = this->read_class_data_.begin();
       data_it != this->read_class_data_.end();
       data_it++) {
       	if (allele_repeats.find(*data_it) == allele_repeats.end()){
       		allele_repeats[*data_it] = 1;
       	}
       	else{
       		allele_repeats[*data_it]++;
       	}
  	}

  	for (map<int32_t, int32_t>::iterator it = allele_repeats.begin(); it != allele_repeats.end(); it++){
  		if (it->second >= 2){
		    (*alleles).push_back(it->first);
  			// cerr << it->first << "\t" << it -> second << endl;
  		}
  	}
	return true;	//TODO add false
}
