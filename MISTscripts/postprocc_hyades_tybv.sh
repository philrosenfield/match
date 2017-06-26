#!/bin/bash

clustn=$MATCH_ODYOUT
jobid=$(ls -t $MATCH_ODYOUT/Hyades/output/hyades_debruTIGS_TychoBV/ | head -1)
jobid=${jobid##*M}
./plotruns.py Hyades hyades_debruTIGS_TychoBV.phot $jobid '$jobid.csv' Tycho_B Tycho_V
