#!/bin/bash
scriptDir=$(cd $(dirname -- $0); pwd)
python ${scriptDir}/Lcfit.py /sharefs/hbkg/data/SCAN/luoqi/HE_GPS_SML/29506101/config_he.xml
