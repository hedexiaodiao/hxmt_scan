import numpy as np
from astropy.io import fits as pf
import sys,os
from glob import glob as glob
sys.path.append("/hxmt/work/USERS/nangyi/tools")
from Fitstogether2 import fitstogetherMulti

listall = []
with open('041crab.all','r')as f:
    for i in f:
        listall.append(i[:13])

attpath = '/sharefs/hbkg/data/SCAN/nangyi/HE_data/att/'
lcpath = '/sharefs/hbkg/data/SCAN/nangyi/HE_data/lc_1811/netlc/'
outpath= '/sharefs/hbkg/data/SCAN/PSF/HE/all_2020/'
attall = [attpath+i+'shortAtt.fits' for i in listall]
lc0all = [lcpath+i+'he_netlc0.fits' for i in listall]
lc1all = [lcpath+i+'he_netlc1.fits' for i in listall]
lc2all = [lcpath+i+'he_netlc5.fits' for i in listall]

fitstogetherMulti(lc0all,outpath+'he_netlc0.fits')
fitstogetherMulti(lc1all,outpath+'he_netlc1.fits')
fitstogetherMulti(lc2all,outpath+'he_netlc5.fits')
fitstogetherMulti(attall,outpath+'crtAtt.fits')


