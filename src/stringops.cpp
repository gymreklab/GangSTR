// Modified from https://github.com/tfwillems/HipSTR/blob/3a7c6f043d822a74df3ad42be00afc35b91f12fa/src/stringops.cpp
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "stringops.h"

void split_by_delim(const std::string &s, char delim, 
		    std::vector<std::string>& substrings){
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim))
    substrings.push_back(item);
}

std::string uppercase(std::string str){
  std::stringstream res;
  for (size_t i = 0; i < str.size(); i++)
    res << static_cast<char>(toupper(str[i]));
  return res.str();
}

std::string lowercase(std::string str){
  std::stringstream res;
  for (size_t i = 0; i < str.size(); i++)
    res << static_cast<char>(tolower(str[i]));
  return res.str();
}

bool string_starts_with(std::string&s, std::string prefix){
  if (s.size() < prefix.size())
    return false;
  return s.substr(0, prefix.size()).compare(prefix) == 0;
}

bool string_ends_with(std::string& s, std::string suffix){
  if (s.size() < suffix.size())
    return false;
  return s.substr(s.size()-suffix.size(), suffix.size()).compare(suffix) == 0;
}

bool orderByLengthAndSequence(const std::string& s1, const std::string s2){
  if (s1.size() != s2.size())
    return s1.size() < s2.size();
  return s1.compare(s2) < 0;
}


int length_suffix_match(std::string& s1, std::string& s2){
  std::string::const_reverse_iterator iter_1 = s1.rbegin();
  std::string::const_reverse_iterator iter_2 = s2.rbegin();
  int num_matches = 0;
  for (; iter_1 != s1.rend() && iter_2 != s2.rend(); ++iter_1, ++iter_2){
    if (*iter_1 != *iter_2)
      break;
    num_matches++;
  }
  return num_matches;
}

std::string reverse_complement(std::string nucs) {
  std::string rev;
  size_t size = nucs.size();
  rev.resize(size);
  for (size_t i = 0; i < size; i++) {
    rev.replace(size-i-1, 1, 1, complement(nucs[i]));
  }
  return rev;
}

char complement(const char nucleotide) {
  switch (nucleotide) {
  case 'A':
  case 'a':
    return 't';
  case 'T':
  case 't':
    return 'a';
  case 'G':
  case 'g':
    return 'c';
  case 'C':
  case 'c':
    return 'g';
  default:
    return 'N';
  }
  return 'N';
}
