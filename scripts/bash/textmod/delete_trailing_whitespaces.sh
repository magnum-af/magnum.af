#!/bin/bash
# enforce exactly one ws after comma (excluding comma followed by newline)
cd $( dirname "${BASH_SOURCE[0]}" ) # cd script's directory (makes it executable form everywhere and avoids unwanted behaviour)

find ../../../../magnum.af/ -name '*.hpp' -type f -exec sed -i -e 's/[[:space:]]\+$//' -- {} +
find ../../../../magnum.af/ -name '*.cpp' -type f -exec sed -i -e 's/[[:space:]]\+$//' -- {} +
find ../../../../magnum.af/ -name '*.py'  -type f -exec sed -i -e 's/[[:space:]]\+$//' -- {} +
find ../../../../magnum.af/ -name '*.pyx' -type f -exec sed -i -e 's/[[:space:]]\+$//' -- {} +
