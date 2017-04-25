#ifndef SRC_COMMON_H__
#define SRC_COMMON_H__

#include <string>
#include <vector>

// Print msg, exit if error
enum MSGTYPE {
  ERROR = 0,
  WARNING,
  DEBUG,
  PROGRESS
};
void PrintMessageDieOnError(const std::string& msg,
                            MSGTYPE msgtype);

// String split
void split_by_delim(const std::string &s, char delim, 
		    std::vector<std::string>& substrings);

#endif  // SRC_COMMON_H__
