#!/bin/bash +e

mode="%p" # test run printing to stdout
# mode=x # apply changes

for file in {*.cpp,*/*.cpp,*.hpp,*/*.hpp}; do
    echo $file
    printf '\n%s\n' '0?#include?a' 'namespace magnumafcpp{' . "$mode" | ex $file
    [ "$mode" == "x" ] && echo "}// namespace magnumafcpp" >> $file
done
