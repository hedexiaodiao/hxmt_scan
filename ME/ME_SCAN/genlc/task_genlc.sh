#!/bin/bash
scriptDir=$(cd $(dirname -- $0); pwd)
source /sharefs/hbkg/user/luoqi/home/scan_MEgenlc_env.sh
python ${scriptDir}/timing_run.py ${1}
