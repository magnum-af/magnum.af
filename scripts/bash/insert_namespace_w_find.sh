#!/bin/bash

mode="%p" # test run printing to stdout
#mode=x # apply changes

files=$(find . -name "*.cpp")
for file in {$files}; do
    echo $file
    printf '\n%s\n' '0?#include?a' 'using namespace magnumaf;' . "$mode" | ex $file
done
