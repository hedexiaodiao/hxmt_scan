import numpy as np
from astropy.io import fits as pf
import sys,os
from glob import glob as glob
###sys.path.append("/hxmt/work/USERS/nangyi/tools")
from Fitstogether2 import fitstogetherMulti

listall = []
with open('041crab.all','r')as f:
    for i in f:
        listall.append(i[:13])

attpath = '/sharefs/hbkg/data/SCAN/luoqi/HE/Org/Att/'
lcpath = '/sharefs/hbkg/data/SCAN/luoqi/HE/Net/'
outpath= '/sharefs/hbkg/data/SCAN/luoqi/HE/PSF/'###'/sharefs/hbkg/data/SCAN/PSF/HE/all_2020/'
attall = [attpath+i+'_shortAtt.fits' for i in listall]
lc0all = [lcpath+i+'/'+i+'_he_netlc_b0.fits' for i in listall]
lc1all = [lcpath+i+'/'+i+'_he_netlc_b1.fits' for i in listall]
lc2all = [lcpath+i+'/'+i+'_he_netlc_b2.fits' for i in listall]

fitstogetherMulti(lc0all,outpath+'he_netlc_b0.fits')
fitstogetherMulti(lc1all,outpath+'he_netlc_b1.fits')
fitstogetherMulti(lc2all,outpath+'he_netlc_b2.fits')
fitstogetherMulti(attall,outpath+'crtAtt.fits')


