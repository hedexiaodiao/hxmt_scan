#!/bin/bash
scriptDir=$(cd $(dirname -- $0); pwd)
hep_sub ${scriptDir}/task_lcfit.sh -argu ${1} -g hxmt -mem 8000 -o /sharefs/hbkg/data/SCAN/HE/joboutput -e /sharefs/hbkg/data/SCAN/HE/joboutput
