#!/bin/bash +e

mode="%p" # test run printing to stdout
#mode=x # apply changes

for file in *.hpp; do
    echo $file
    printf '\n%s\n' '0?#include?a' 'namespace magnumaf{' . "$mode" | ex $file
    [ "$mode" == "x" ] && echo "}// namespace magnumaf" >> $file
done

for file in */*.hpp; do
    echo $file
    printf '\n%s\n' '0?#include?a' 'namespace magnumaf{' . "$mode" | ex $file
    [ "$mode" == "x" ] && echo "}// namespace magnumaf" >> $file
done

for file in *.cpp; do
    echo $file
    #printf '\n%s\n' '0?#include?a' 'using namespace magnumaf;' . "$mode" | ex $file
    printf '\n%s\n' '0?#include?a' 'namespace magnumaf{' . "$mode" | ex $file
    [ "$mode" == "x" ] && echo "}// namespace magnumaf" >> $file
done

for file in */*.cpp; do
    echo $file
    #printf '\n%s\n' '0?#include?a' 'using namespace magnumaf;' . "$mode" | ex $file
    printf '\n%s\n' '0?#include?a' 'namespace magnumaf{' . "$mode" | ex $file
    [ "$mode" == "x" ] && echo "}// namespace magnumaf" >> $file
done
