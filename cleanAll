#!/bin/bash

# ===
# Script for cleaning all of the tutorials.
# ===

for tutorialDir in */ ; do
    if [[ $tutorialDir = *OFtutorial* ]]; then
        echo "Cleaning:" $tutorialDir

        cd $tutorialDir
        ./Allwclean >/dev/null 2>&1
        cd testCase
        ./Allclean >/dev/null 2>&1
        cd ../..
    fi
done
