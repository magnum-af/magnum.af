#ifndef MISC_H
#define MISC_H
#include <sys/stat.h>
#include <unistd.h>
#include <string>

unsigned long long GetDirSize(std::string filepath);
inline bool exists (const std::string& absolute_filepath);

///> Returns true if file "name" exists, false otherwise
inline bool exists (const std::string& absolute_filepath) {
    struct stat buffer;   
    return (stat (absolute_filepath.c_str(), &buffer) == 0); 
}

///> Colorizing a std::string for std::cout. Red is for warnings.
std::string red(const std::string str);
///> Colorizing a std::string for std::cout. Bold red is for errors.
std::string bold_red(const std::string str);
///> Colorizing a std::string for std::cout. Green is for regular infos.
std::string bold_green(const std::string str);
///> Colorizing a std::string for std::cout. Bold green is for sucess infos.
std::string green(const std::string str);
#endif
