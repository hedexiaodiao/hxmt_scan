#!/bin/bash
source /sharefs/hbkg/user/luoqi/home/mypython
scriptDir=$(cd $(dirname -- $0); pwd)
python ${scriptDir}/genlc_module.py ${1}
