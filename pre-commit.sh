#!/bin/sh
git stash -q --keep-index
./install
./run_short_tests
RESULT=$?
git stash pop -q
[ $RESULT -ne 0 ] && exit 1
exit 0