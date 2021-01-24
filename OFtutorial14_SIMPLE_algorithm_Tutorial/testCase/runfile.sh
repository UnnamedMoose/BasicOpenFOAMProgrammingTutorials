#!/bin/bash

rm -f log.*
foamListTimes -rm
foamCleanPolyMesh
blockMesh 2>&1 | tee log.blockmesh
simpleFoamTutorial 2>&1 | tee log.simulation
