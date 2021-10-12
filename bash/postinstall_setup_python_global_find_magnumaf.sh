#!/bin/bash

# This script adds a .pth file such that python finds the installed library mangumaf.so for every user. The .pth file contains the path to the magnumaf.so library

echo '/usr/local/lib' | sudo tee /usr/lib/python3/dist-packages/find_magnumaf.pth
