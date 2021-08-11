#!/bin/bash -e
# assumes dir ../ is project root dir

# call this scripts directory
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# cd into project root dir
cd ../

echo "Running formatting from dir $PWD"

set -o xtrace

# format C++
find . -iname '*.hpp' -o -iname '*.cpp' | xargs clang-format -i --style=file

# format CMake
find . -iname 'CMakeLists.txt' | xargs cmake-format -i --line-width=120
