#!/bin/bash
# $1 old pattern
# $2 new string
git grep -l '$1' | xargs sed -i 's/$1/$2/g'
