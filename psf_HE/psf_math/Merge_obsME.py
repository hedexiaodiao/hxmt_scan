import numpy as np
from astropy.io import fits as pf
import sys,os
from glob import glob as glob
###sys.path.append("/hxmt/work/USERS/nangyi/tools")
from Fitstogether2 import fitstogetherMulti

lc0all = ['A01_me_b0.fits','A02_me_b0.fits']
lc1all = ['A01_me_b1.fits','A02_me_b1.fits']
lc2all = ['A01_me_b2.fits','A02_me_b2.fits']
attall = ['A01_me_att.fits','A02_me_att.fits']
###print(attall)
fitstogetherMulti(lc0all,'./A01A02_me_b0.fits')
fitstogetherMulti(lc1all,'./A01A02_me_b1.fits')
fitstogetherMulti(lc2all,'./A01A02_me_b2.fits')
fitstogetherMulti(attall,'./A01A02_me_att.fits')


