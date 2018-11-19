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

#ifndef SRC_VCF_WRITER_H
#define SRC_VCF_WRITER_H__

#include <iostream>
#include <fstream>

#include "src/locus.h"

using namespace std;

class VCFWriter {
 public:
  VCFWriter(const std::string& _vcffile, const std::string& full_command, const std::string& sample_name);
  void WriteRecord(const Locus& locus);
  virtual ~VCFWriter();
 private:
  ofstream writer_;
};

#endif  // SRC_VCF_WRITER_H__
