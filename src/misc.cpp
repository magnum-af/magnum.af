#include <string>
#include <iostream>

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
