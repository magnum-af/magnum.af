#ifndef MISC_H
#define MISC_H
#include <sys/stat.h>
#include <unistd.h>
#include <string>

unsigned long long GetDirSize(std::string filepath);
inline bool exists (const std::string& absolute_filepath);

// Returns true if file "name" exists, false otherwise
inline bool exists (const std::string& absolute_filepath) {
    struct stat buffer;   
    return (stat (absolute_filepath.c_str(), &buffer) == 0); 
}

#endif
