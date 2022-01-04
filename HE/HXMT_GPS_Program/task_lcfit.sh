#!/bin/bash
source /sharefs/hbkg/user/luoqi/home/mypython
scriptDir=$(cd $(dirname -- $0); pwd)
python ${scriptDir}/Lcfit.py ${1}
