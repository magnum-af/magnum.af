#!/bin/bash +e

mode="%p" # test run printing to stdout
#mode=x # apply changes

for file in *.hpp; do
    echo $file
    printf '\n%s\n' '0?#include?a' 'namespace magnumafcpp{' . "$mode" | ex $file
    [ "$mode" == "x" ] && echo "}// namespace magnumafcpp" >> $file
done

for file in */*.hpp; do
    echo $file
    printf '\n%s\n' '0?#include?a' 'namespace magnumafcpp{' . "$mode" | ex $file
    [ "$mode" == "x" ] && echo "}// namespace magnumafcpp" >> $file
done

for file in *.cpp; do
    echo $file
    #printf '\n%s\n' '0?#include?a' 'using namespace magnumafcpp;' . "$mode" | ex $file
    printf '\n%s\n' '0?#include?a' 'namespace magnumafcpp{' . "$mode" | ex $file
    [ "$mode" == "x" ] && echo "}// namespace magnumafcpp" >> $file
done

for file in */*.cpp; do
    echo $file
    #printf '\n%s\n' '0?#include?a' 'using namespace magnumafcpp;' . "$mode" | ex $file
    printf '\n%s\n' '0?#include?a' 'namespace magnumafcpp{' . "$mode" | ex $file
    [ "$mode" == "x" ] && echo "}// namespace magnumafcpp" >> $file
done
