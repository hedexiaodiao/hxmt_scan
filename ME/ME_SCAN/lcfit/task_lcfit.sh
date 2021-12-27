#!/bin/bash
source /sharefs/hbkg/user/luoqi/home/scan_MEfit_env.sh
scriptDir=$(cd $(dirname -- $0); pwd)
python ${scriptDir}/test_all_typeSub.py ${1}
