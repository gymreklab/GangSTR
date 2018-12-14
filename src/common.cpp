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

#include <err.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>

#include "src/common.h"

using namespace std;

void PrintMessageDieOnError(const string& msg, MSGTYPE msgtype) {
  string typestring = "";
  switch (msgtype) {
  case M_ERROR:
    typestring = "ERROR: ";
    break;
  case M_WARNING:
    typestring = "WARNING: ";
    break;
  case M_PROGRESS:
    typestring = "ProgressMeter: ";
    break;
  case M_DEBUG:
    typestring = "DEBUG: ";
    break;
  default:
    errx(1, "Invalid message type. This should never happen");
  }
  stringstream ss;
  ss  << "[GangSTR"
      << "-" << _GIT_VERSION << "] " << typestring << msg << endl;
  cerr << ss.str();

  if (msgtype == M_ERROR) {
    exit(1);
  }
}

float GetGC(const std::string& seq) {
  int total_bases = 0;
  int gc = 0;
  for (int i=0; i<seq.size(); i++) {
    if (seq[i] != 'c' && seq[i] != 'C' &&
	seq[i] != 'g' && seq[i] != 'G' &&
	seq[i] != 't' && seq[i] != 'T' &&
	seq[i] != 'a' && seq[i] != 'A') {
      return -1;
    }
    total_bases += 1;
    if (seq[i] == 'c' || seq[i] =='C' ||
	seq[i] == 'g' || seq[i] == 'G') {
      gc += 1;
    }
  }
  return float(gc)/float(total_bases);
}
