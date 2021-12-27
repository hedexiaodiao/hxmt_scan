#!/bin/bash
scriptDir=$(cd $(dirname -- $0); pwd)
hep_sub ${scriptDir}/task_genlc.sh -argu P010129506101 -g hxmt -mem 8000 -o /sharefs/hbkg/data/SCAN/luoqi/joboutput -e /sharefs/hbkg/data/SCAN/luoqi/joboutput
