#!/bin/bash
scriptDir=$(cd $(dirname -- $0); pwd)
python ${scriptDir}/timing_run.py ${1}
