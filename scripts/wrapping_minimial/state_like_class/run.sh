#!/bin/bash

python setup.py build_ext -i
LD_LIBRARY_PATH=/usr/local/arrayfire/lib python run.py

