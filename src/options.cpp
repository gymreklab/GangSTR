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

using namespace std;

Options::Options() {
  // initialize options
  reffa = "";
  regionsfile = "";
  outprefix = "";
  dist_mean = 400;
  dist_sdev = 50;
  stutter_up = 0.01;
  stutter_down = 0.02;
  stutter_p = 0.95;
  flanklen = 2000;
  realignment_flanklen = 100;
  frr_weight = 0.8;
  enclosing_weight = 1.0;
  spanning_weight = 1.0;
  verbose = false;
}

Options::~Options() {}
