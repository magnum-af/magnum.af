#!/bin/bash
# this scripts tests wether directory in $1 exists

# determining write directory
if [ -n "$1" ];then
    if [ -e "$1" ];then
        echo "Error: Write Directory exists!" $1
	exit 1
    else
        mkdir --parents $1
        echo "writing in new directory " $1
        exit 0
    fi
fi
