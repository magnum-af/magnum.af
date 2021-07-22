#!/bin/bash +e

mode="%p" # test run printing to stdout
# mode=x # apply changes

for file in {*.cpp,*/*.cpp,*.hpp,*/*.hpp}; do
    echo $file
    printf '\n%s\n' '0?#include?a' 'namespace magnumaf{' . "$mode" | ex $file
    [ "$mode" == "x" ] && echo "}// namespace magnumaf" >> $file
done
