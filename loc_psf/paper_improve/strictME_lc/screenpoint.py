from astropy.io import fits as pf
import numpy as np
import glob
import sys,os
sys.path.append("/sharefs/hbkg/user/nangyi/tools")
from Fitstogether2 import *

def distance(ratpp,dectpp,rapnt,decpnt):
    ytp = (np.sin(decpnt*np.pi/360.0 - dectpp*np.pi/360.0))**2 + np.cos(decpnt*np.pi/180.0)*np.cos(dectpp*np.pi/180.0)*(np.sin(ratpp*np.pi/360.0 - rapnt*np.pi/360.0))**2
    ytp = np.arccos(1-2*ytp)*180./np.pi
    return ytp

def Uniform(arrylist):
    '''
    Input list (N dimensions),return the intersection of N dimension data
    '''
    if arrylist.count(arrylist[0]) == len(arrylist):
        return arrylist[0]
    times = sum(arrylist,[])
    lenth = len(arrylist)
    timeset = list(set(times))
    unitime = [i for i in timeset if times.count(i)==lenth]
    unitime.sort()
    return unitime

def genfits(lis,outname):
    col1 = pf.Column(name = 'Time', format='D',array = lis[0])
    col2 = pf.Column(name = 'Counts', format='D',array = lis[1])
    col3 = pf.Column(name = 'Stat_err', format='D',array = lis[2])
    cols = pf.ColDefs([col1,col2,col3])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    tbhdu.writeto(outname,clobber=True)

lepath = '/sharefs/hbkg/data/SCAN/nangyi/LE_data/lc/'
mepath = '/sharefs/hbkg/data/SCAN/nangyi/ME_data/lc/'
hepath = '/sharefs/hbkg/data/SCAN/nangyi/HE_point/lc/'

ledatas = glob.glob(lepath+'P0101299*lebox0.fits')
medatas = glob.glob(mepath+'P0101299*mebox0.fits')
hedatas = glob.glob(hepath+'P0101299*hebox0.fits')

ledatas=[i.split('/')[-1][:13] for i in ledatas]
medatas=[i.split('/')[-1][:13] for i in medatas]
hedatas=[i.split('/')[-1][:13] for i in hedatas]

unidatas = Uniform([ledatas,medatas,hedatas])
outpath = '/sharefs/hbkg/data/SCAN/nangyi/pointcrab3/'
for i in range(3):
    os.system('rm %s; rm %s; rm %s'%(outpath+'lepoint%s_all.fits'%i,outpath+'mepoint%s_all.fits'%i,outpath+'hepoint%s_all.fits'%i))
    ledatas = [lepath+j+'lebox%s.fits'%i for j in unidatas]
    fitstogetherMulti(ledatas,outpath+'lepoint%s_all.fits'%i)
    medatas = [mepath+j+'mebox%s.fits'%i for j in unidatas]
    fitstogetherMulti(medatas,outpath+'mepoint%s_all.fits'%i)
    hedatas = [hepath+j+'hebox%s.fits'%i for j in unidatas]
    fitstogetherMulti(hedatas,outpath+'hepoint%s_all.fits'%i)

atts = ['/sharefs/hbkg/data/SCAN/nangyi/HE_point/att/'+j+'crtatt' for j in unidatas]
os.system('rm %s'%(outpath+'pointatt.fits'))
fitstogetherMulti(atts,outpath+'pointatt.fits')

outpath2 = '/sharefs/hbkg/data/SCAN/nangyi/pointcrab3/'

hd=pf.open(outpath+'hepoint0_all.fits')
hct = hd[1].data.field(0)
hd.close();del hd;
hd=pf.open(outpath+'mepoint0_all.fits')
mct = hd[1].data.field(0)
hd.close();del hd;
hd=pf.open(outpath+'lepoint0_all.fits')
lct = hd[1].data.field(0)
hd.close();del hd;
hd=pf.open(outpath+'pointatt.fits')
attt = hd[1].data.field(0)
hd.close();del hd;
attm=attt[np.in1d(attt,lct)]
unitm = Uniform([hct.tolist(),mct.tolist(),lct.tolist(),attm.tolist()])
#unitm=unitm[100:30000]
#hd=pf.open(outpath2+'shortpointatt.fits')
#unitm = hd[1].data.field(0)
print len(unitm)
#hd.close();del hd
for i in range(3):
    hd=pf.open(outpath+'hepoint%s_all.fits'%i)
    tm = hd[1].data.field(0)
    yp=np.in1d(tm,unitm)
    hd[1].data=hd[1].data[yp]
    hd.writeto(outpath2+'hepoint%s.fits'%i,overwrite='True')
    hd.close();del hd
    hd=pf.open(outpath+'mepoint%s_all.fits'%i)
    tm = hd[1].data.field(0)
    yp=np.in1d(tm,unitm)
    hd[1].data=hd[1].data[yp]
    hd.writeto(outpath2+'mepoint%s.fits'%i,overwrite='True')
    hd.close();del hd
    hd=pf.open(outpath+'lepoint%s_all.fits'%i)
    tm = hd[1].data.field(0)
    yp=np.in1d(tm,unitm)
    hd[1].data=hd[1].data[yp]
    hd.writeto(outpath2+'lepoint%s.fits'%i,overwrite='True')
    hd.close();del hd

atthd=pf.open(outpath+'pointatt.fits')
attm = atthd[1].data.field(0)
yep = np.in1d(attm,unitm)
for i in range(1,len(atthd),1):
    atthd[i].data = atthd[i].data[yep]

atthd.writeto(outpath2+'shortpointatt.fits',overwrite=True)
atthd.close()

os.system('rm %spointinfo*.txt; rm %s?epoint?.???'%(outpath,outpath))
