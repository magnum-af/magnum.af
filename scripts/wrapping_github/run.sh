#!/bin/bash
rm interface.cpp && python setup.py build_ext -i && python run.py
