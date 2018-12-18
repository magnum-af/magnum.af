#!/bin/bash
if [ -f interface.cpp ];then
    rm interface.cpp
fi
if [ -d ./build ];then
    rm -r ./build
fi
python setup.py build_ext -i 
python run.py

# cleanup
if [ -f interface.cpp ];then
    rm interface.cpp
fi
if [ -d ./build ];then
    rm -r ./build
fi
if [ -f ./wrap.so ];then
    rm ./wrap.so
fi
