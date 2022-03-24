from astropy.io import fits as pf
import os
import numpy as np
import glob
path = '/sharefs/hbkg/data/SCAN/nangyi/ME_data/gti/'
outpath = '/sharefs/hbkg/data/SCAN/locME/gti/'
lts = glob.glob(path+'*mebd.fits')#+glob.glob('/sharefs/hbkg/user/saina/data294/P0101294/*/bd_P*.fits')

badnum=set({})
a=np.array([])
for i in lts:
    hd=pf.open(i)
    tmp = hd[1].data.field(0)
    print (tmp<576).sum()
    badnum = badnum.union(tmp)
    hd.close();del hd


###hd=pf.open('/sharefs/hbkg/data/SCAN/nangyi/ME_data/gti/P010132800302mebd.fits')
hd=pf.open('/sharefs/hbkg/data/SCAN/nangyi/ME_data/gti/P010132800803mebd.fits')
prihdr=hd[1].header

prihdu = pf.PrimaryHDU(header=prihdr)

lth = len(badnum)
col1 = pf.Column(name = 'DetID', format='I',array = np.sort(list(badnum)))
col2 = pf.Column(name = 'TIMERANGE', format='20A',array = np.array(['0']*lth))
col3 = pf.Column(name = 'TIMERANGE2', format='20A',array = np.array(['INDEF']*lth))
col4 = pf.Column(name = 'TYPE', format='20A',array = np.array(['Bad']*lth))
col5 = pf.Column(name = 'STATUS', format='B',array = np.array([0]*lth))
cols = pf.ColDefs([col1,col2,col3,col4,col5])
tbhdu = pf.BinTableHDU.from_columns(cols)
pf.conf.extension_name_case_sensitive=True
tbhdu.name='detectorStatu'
tbhdu.header = prihdr

tbhdu.writeto(outpath + "mebadall.fits",overwrite=True)


