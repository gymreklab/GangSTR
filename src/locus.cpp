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

#include "src/locus.h"

using namespace std;

Locus::Locus() {
  chrom = "";
  start = -1;
  end = -1;
  period = -1;

  insert_size_mean = -1.0;
  insert_size_stddev = -1.0;
  allele1 = -1;
  allele2 = -1;
  lob1 = -1;
  hib1 = -1;
  lob2 = -1;
  hib2 = -1;
  min_neg_lik = 0;
}

Locus::~Locus() {}
