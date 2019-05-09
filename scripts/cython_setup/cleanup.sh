#!/bin/bash
cd "$(dirname "${BASH_SOURCE[0]}")"
[[ -d "build/" ]] && rm -r build/ && echo "removed build/"
deliting="magnumaf.cpp magnum_af_cython_setup.so magnum_af_cython_setup.cpython-36m-x86_64-linux-gnu.so magnumaf_decl.pxd  magnumaf.pyx"
for file in $deliting; do
    [[ -f "$file" ]] && rm "$file" && echo "removed $file"
done
