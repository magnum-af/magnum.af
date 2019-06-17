#ifndef MISC_H
#define MISC_H
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <cstring>
#include <pwd.h>
#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

unsigned long long GetDirSize(std::string filepath);
///> Returns true if file "name" exists, false otherwise
inline bool exists (const std::string& absolute_filepath) {
    struct stat buffer;
    return (stat (absolute_filepath.c_str(), &buffer) == 0);
}

inline bool createdir(const std::string& absolute_filepath);
std::string setup_magafdir();

//void directoryManager(char *dir, int maxNumberOfFiles);
void remove_oldest_files_until_size(const char *dir, unsigned long long  maxNumberOfBytes, bool verbose = true);

///> Colorizing a std::string for std::cout. Red is for warnings.
inline std::string red(const std::string str){return "\033[;31m"+str+"\033[0m";}
///> Colorizing a std::string for std::cout. Bold red is for errors.
inline std::string bold_red(const std::string str){return "\033[1;31m"+str+"\033[0m";}
///> Colorizing a std::string for std::cout. Green is for regular infos.
inline std::string bold_green(const std::string str){return "\033[1;32m"+str+"\033[0m";}
///> Colorizing a std::string for std::cout. Bold green is for sucess infos.
inline std::string green(const std::string str){return "\033[;32m"+str+"\033[0m";}

inline const char* Info(void){ return "\33[0;32mInfo:\33[0m";}
inline const char* Warning(void){ return "\33[1;31mWarning:\33[0m";}

#endif
