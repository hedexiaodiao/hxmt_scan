#!/bin/bash
scriptDir=$(cd $(dirname -- $0); pwd)
hep_sub ${scriptDir}/task_allpro_lc_only.sh -argu ${1} -g hxmt -mem 9000 -o /sharefs/hbkg/data/SCAN/HE/joboutput -e /sharefs/hbkg/data/SCAN/HE/joboutput
