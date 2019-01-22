#!/bin/bash
# this scripts tests wether directory in $1 exists
# $2 ... verbose flag

# determining write directory
if [ -n "$1" ];then
    if [ -e "$1" ];then
        echo "Error: Write Directory '$1' exists! Use the '-f' flag to overwrite"
	exit 1
    else
        if [ -d "$1" ];then
            mkdir --parents $1
        else
            mkdir --parents $1/
        fi
        if [ "$2" == "true" ];then
            echo "writing in new directory " $1
        fi
        exit 0
    fi
fi
