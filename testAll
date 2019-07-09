#!/bin/bash

# ===
# Script for automatically compiling, running and testing all of the tutorials.
# ===

# run a command and redirect outputs to /dev/null, only print PASS or FAIL
# depending on whether errors occured or not
run () {
    # run the command and get rid of the outputs
    $1 >/dev/null 2>&1
    # check return code, 0 indicates success
    if [ $? -eq 0 ]; then
        echo "    PASS: "$1
    else
        echo "    FAIL: "$1
    fi
}

# go over all tutorial directories
for tutorialDir in */ ; do
    if [[ $tutorialDir = *OFtutorial* ]]; then
        echo "Checking:" $tutorialDir

        # navigate to the tutorial
        cd $tutorialDir

        # no need to test cleaning the code
        ./Allwclean >/dev/null 2>&1

        # test building and running
        run ./Allwmake
        cd testCase
        run ./Allrun

        # no need to test cleaning the test case
        ./Allclean >/dev/null 2>&1

        # go back to main directory
        cd ../..
    fi
done
