#!/bin/bash -e
# running all cpp integration test binaries
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # call this scripts directory
for file in bin/*; do $file; done
