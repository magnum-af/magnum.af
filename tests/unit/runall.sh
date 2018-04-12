#!/bin/bash
PYTHONPATH=../../build/src/ python atomistic_anisotropy.py
PYTHONPATH=../../build/src/ python atomistic_dipole_dipole.py
PYTHONPATH=../../build/src/ python atomistic_exchange.py
PYTHONPATH=../../build/src/ python atomistic_dmi.py

#Planned:
#tests for all interactions
#tests for integration methods during refactoring
#"integration tests" for sp4, maybe sp5
