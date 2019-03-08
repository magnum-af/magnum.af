#include "misc.hpp"
#include <iostream>
#include <string>

unsigned long long GetDirSize(std::string filepath)
// Retuns size of directory in bytes
// filepath is expected to be an absolute file path as the first word is detected by the first occurance of "/"
// http://5eonline.com/en/4-easy-ways-to-run-external-programs-in-cc/
{
    FILE *stream;
    char string[ 256 ] = { 0 };
    stream = popen(("du --max-depth=0 "+filepath).c_str(), "r" );
    while( NULL != fgets( string, sizeof( string ) - 1, stream ) )
    {
       //printf( "(%s)\n", string );
    }
    std::string charstring(string);
    std::string dir_size_string=charstring.substr(0,charstring.find("/")-1);
    return std::stoull (dir_size_string,0,0);

}

inline bool createdir(const std::string& absolute_filepath){
// Creating a directory
    if (mkdir(absolute_filepath.c_str(), 0777) == -1){
        std::cerr << "Error in createdir for " << absolute_filepath<< " : " << std::strerror(errno) << std::endl;
        return false;
    }

    else{
        //printf("Directory ~/.magnum_af/ created.");
        std::cout << "Directory"+absolute_filepath+"created" << std::endl;
        return true;
    }

}

std::string setup_magafdir(){
    const char *homedir;
    if ((homedir = getenv("HOME")) == NULL) {
        homedir = getpwuid(getuid())->pw_dir;
    }
    std::string magafdir(homedir);
    magafdir = magafdir +"/.magnum.af.cache/";
    if (exists(magafdir) == false){
       createdir(magafdir); 
    }
    return magafdir;
}

std::string red(const std::string str){
    return "\033[;31m"+str+"\033[0m";
}

std::string bold_red(const std::string str){
    return "\033[1;31m"+str+"\033[0m";
}

std::string bold_green(const std::string str){
    return "\033[1;32m"+str+"\033[0m";
}

std::string green(const std::string str){
    return "\033[;32m"+str+"\033[0m";
}


