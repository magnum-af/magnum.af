#!/bin/bash
hits=$(git grep -l LLGIntegrator -- '*.py' LLGIntegrator)
for file in $hits; do
    echo $file
    vim -c ":%s/LLGIntegrator(\(.*\))/LLGIntegrator([\1])/g" $file -c x
done
