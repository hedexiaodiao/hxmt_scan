#!/bin/bash
source /sharefs/hbkg/user/luoqi/home/mypython
scriptDir=$(cd $(dirname -- $0); pwd)
python ${scriptDir}/funcgenlc_module_lc_only.py ${1}
