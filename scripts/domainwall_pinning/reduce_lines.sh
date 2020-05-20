#!/bin/bash -ex
# $1 filename
# $2 optional: every line only

filename="${1-"m.dat"}"
every_line_only="${2-1000}"
awk '(NR%'$every_line_only'==1)' "$filename" > "$filename".reduced
tail -n 1 "$filename" >> "$filename".reduced
