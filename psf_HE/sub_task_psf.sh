#!/bin/bash
scriptDir=$(cd $(dirname -- $0); pwd)
hep_sub ${scriptDir}/task_psf.sh -argu ${1} -g hxmt -mem 8000 -o /sharefs/hbkg/data/SCAN/luoqi/joboutput -e /sharefs/hbkg/data/SCAN/luoqi/joboutput