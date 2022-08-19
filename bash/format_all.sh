#!/bin/bash -e
# assumes dir ../ is project root dir

# call this scripts directory
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# cd into project root dir
cd ../

echo "Running formatting from dir $PWD"

set -o xtrace

# format C++
find magnumaf/ -iname '*.hpp' -o -iname '*.cpp' \
    | xargs -P 0 -n 8 clang-format -i --style=file --verbose

# format CMake
find . -iname 'CMakeLists.txt' \
    | xargs -P 0 cmake-format -i --line-width=120

# format python
autopep8 python/ --recursive --in-place --verbose
