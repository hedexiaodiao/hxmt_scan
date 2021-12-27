#!/bin/bash
scriptDir=$(cd $(dirname -- $0); pwd)
hep_sub ${scriptDir}/task_genlc.sh -argu ${1} -g hxmt -mem 9000 -o /sharefs/hbkg/data/SCAN/HE/joboutput -e /sharefs/hbkg/data/SCAN/HE/joboutput
