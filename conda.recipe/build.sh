#!/bin/bash

# For debugging
#set -o verbose

# Die on unset variables
set -o nounset
# Die on errors
set -o errexit
# Die if any part of a pipe fails
set -o pipefail

python setup.py install --single-version-externally-managed --record=record.txt
make -C src/ libpyhuddingec.so
cp src/libpyhuddingec.so ${PREFIX}/lib/libpyhuddingec.so