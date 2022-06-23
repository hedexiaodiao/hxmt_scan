import numpy as np
from astropy.io import fits as pf
import sys,os
from glob import glob as glob
###sys.path.append("/hxmt/work/USERS/nangyi/tools")
from Fitstogether2 import fitstogetherMulti


outpath = '.'

listall = []
with open('crablist.txt','r')as f:
    for i in f:
        listall.append(i[:13])

print(listall)


attall = []
lc0all = []
lc1all = []
lc2all = []

for i in range(len(listall)):
    ObsID = listall[i]
    stObsID = ObsID[0:8]

    Org_attpath = '/sharefs/hbkg/data/SCAN/ME/Org/Att'
    Wcutattfile = "%s/%s_cut_Att.fits" % (Org_attpath, ObsID)

    attfile = Wcutattfile
    lc0file = '/sharefs/hbkg/data/SCAN/ME/Net/%s/me_lc_box0_small_cut.fits'%(ObsID)
    lc1file = '/sharefs/hbkg/data/SCAN/ME/Net/%s/me_lc_box1_small_cut.fits'%(ObsID)
    lc2file = '/sharefs/hbkg/data/SCAN/ME/Net/%s/me_lc_box2_small_cut.fits'%(ObsID)
    attflag = os.path.exists(attfile)
    lc0flag = os.path.exists(lc0file)
    lc1flag = os.path.exists(lc1file)
    lc2flag = os.path.exists(lc2file)
    if attflag and lc0flag and lc1flag and lc2flag:
        attall.append(attfile)
        lc0all.append(lc0file)
        lc1all.append(lc1file)
        lc2all.append(lc2file)
    else:
        print("lc file no exist:",ObsID,attflag,lc0flag,lc1flag,lc2flag)

###print(attall)
fitstogetherMulti(lc0all,outpath+'/me_b0.fits')
fitstogetherMulti(lc1all,outpath+'/me_b1.fits')
fitstogetherMulti(lc2all,outpath+'/me_b2.fits')
fitstogetherMulti(attall,outpath+'/me_att.fits')


