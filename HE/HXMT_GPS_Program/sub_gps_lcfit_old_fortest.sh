#!/bin/bash
scriptDir=$(cd $(dirname -- $0); pwd)
hep_sub ${scriptDir}/task_lcfit.sh -argu /sharefs/hbkg/data/SCAN/luoqi/HE_GPS_SML/29506101/config_he.xml -g hxmt -mem 8000 -o /sharefs/hbkg/data/SCAN/luoqi/joboutput -e /sharefs/hbkg/data/SCAN/luoqi/joboutput
