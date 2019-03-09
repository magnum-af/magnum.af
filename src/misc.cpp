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
        //std::cerr << "Error in createdir for " << absolute_filepath<< " : " << std::strerror(errno) << std::endl;
        printf("Error in createdir for %s : %s \n<< ", absolute_filepath.c_str(), std::strerror(errno));
        return false;
    }

    else{
        //printf("Directory ~/.magnum_af/ created.");
        //std::cout << "Directory"+absolute_filepath+"created" << std::endl;
        printf("Directory '%s' created \n", absolute_filepath.c_str());
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

//adapted from https://stackoverflow.com/questions/9642145/is-there-a-way-to-find-the-oldest-file-using-just-the-c
void remove_oldest_files_until_size(const char *dir, unsigned long long  maxNumberOfBytes, bool verbose){
    int maxiter = 0;
    while(GetDirSize(std::string(dir)) >= maxNumberOfBytes && maxiter < 10){

        DIR *dp;
        struct dirent *entry, *oldestFile=NULL;
        struct stat statbuf;
        int numberOfEntries=0;
        time_t t_oldest;

        time(&t_oldest);
        if((dp = opendir(dir)) != NULL) {
             if(chdir(dir)==0){
                 while((entry = readdir(dp)) != NULL) {
                    lstat(entry->d_name, &statbuf);
                    if(strcmp(".",entry->d_name) == 0 || strcmp("..",entry->d_name) == 0)
                       continue;
                    if (maxiter == 0 && verbose) printf("Entry: ~/.magnum.af.cache/%s\t%s", entry->d_name, ctime(&statbuf.st_mtime));
                       numberOfEntries++;
                    if(difftime(statbuf.st_mtime, t_oldest) < 0){
                          t_oldest = statbuf.st_mtime;
                       oldestFile = entry;
                    }
                 }
             }
        }
        if (verbose) printf("Removing oldest file '~/.magnum.af.cache/%s'  %s'\n", oldestFile->d_name, ctime(&t_oldest));
        remove(oldestFile->d_name);
        closedir(dp);
        maxiter++;
    }
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


