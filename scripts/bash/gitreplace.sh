#!/bin/bash
# $1 old pattern
# $2 new string
# Example: gitreplace.sh \\bOldWord\\b NewWord
# Note   : for replacing exact words (and not every pattern) use '\\b' (i.e. '\'-escaped '\b').
set -x
git grep -l "$1" | xargs sed -i "s/$1/$2/g"
