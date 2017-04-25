#ifndef SRC_COMMON_H__
#define SRC_COMMON_H__

// Print msg, exit if error
enum MSGTYPE {
  ERROR = 0,
  WARNING,
  DEBUG,
  PROGRESS
};
void PrintMessageDieOnError(const std::string& msg,
                            MSGTYPE msgtype);

#endif  // SRC_COMMON_H__
