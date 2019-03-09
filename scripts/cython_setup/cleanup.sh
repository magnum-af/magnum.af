#!/bin/bash
cd "$(dirname "${BASH_SOURCE[0]}")"
[[ -d "build/" ]] && rm -r build/ && echo "removed build/"
deliting="magnum_af.cpp magnum_af_cython_setup.so magnum_af_cython_setup.cpython-36m-x86_64-linux-gnu.so magnum_af_decl.pxd  magnum_af.pyx"
for file in $deliting; do
    [[ -f "$file" ]] && rm "$file" && echo "removed $file"
done
