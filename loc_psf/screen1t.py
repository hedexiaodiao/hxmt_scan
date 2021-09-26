from astropy.io import fits as pf
import numpy as np
import glob
import sys,os

#from scipy.signal import find_peaks
sys.path.append("/sharefs/hbkg/user/nangyi/tools")
from Fitstogether2 import *


def FindGap(inptm,gap=1):
    times = np.copy(inptm)
    a=((times[1:]-times[:-1])>gap)
    c=np.nonzero(a)[0]
    potstart = np.r_[0,c+1]
    potstop = np.r_[c,times.shape[0]-1]
    return potstart,potstop

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
outpath = '/sharefs/hbkg/data/SCAN/nangyi/scox1_3/'
outpath = '/sharefs/hbkg/data/SCAN/nangyi/pointcrab1/'

medatas0 = glob.glob(mepath+'P030229000*mebox0.fits') + glob.glob(mepath+'P020204125*mebox0.fits')#glob.glob(mepath+'P0101328*mebox0.fits')
medatas0=sorted([i.split('/')[-1][:13] for i in medatas0])

medatas1 = glob.glob(mepath+'P030229000*mebox0.fits') + glob.glob(mepath+'P020204125*mebox0.fits')#glob.glob(mepath+'P0101328*mebox1.fits')
medatas1=sorted([i.split('/')[-1][:13] for i in medatas1])

medatas2 = glob.glob(mepath+'P030229000*mebox0.fits') + glob.glob(mepath+'P020204125*mebox0.fits')#glob.glob(mepath+'P0101328*mebox2.fits')
medatas2=sorted([i.split('/')[-1][:13] for i in medatas2])

unidatas = sorted([i for i in set(medatas0+medatas1+medatas2)])

os.system('rm %s; rm %s; rm %s'%(outpath+'lepoint*.fits',outpath+'mepoint*.fits',outpath+'hepoint*.fits'))
medatas = [mepath+j+'mebox0.fits' for j in medatas0]
fitstogetherMulti(medatas,outpath+'mepoint0_all.fits')

medatas = [mepath+j+'mebox1.fits' for j in medatas1]
fitstogetherMulti(medatas,outpath+'mepoint1_all.fits')

medatas = [mepath+j+'mebox2.fits' for j in medatas2]
fitstogetherMulti(medatas,outpath+'mepoint2_all.fits')

atts = ['/sharefs/hbkg/data/SCAN/nangyi/ME_data/att/'+j+'crtatt' for j in unidatas]
os.system('rm %s'%(outpath+'pointatt.fits'))
fitstogetherMulti(atts,outpath+'pointatt.fits',noredundant=0)

hdatt = pf.open(outpath+'pointatt.fits')
attm = hdatt[1].data.field(0)
atm_s,order = np.unique(attm,return_index=True)

hdatt[1].data = hdatt[1].data[order]
hdatt[2].data = hdatt[2].data[order]
hdatt[3].data = hdatt[3].data[order]
hdatt[4].data = hdatt[4].data[order]
hdatt.writeto(outpath+'pointatt.fits',overwrite=True)
hdatt.close();del hdatt
os.system('rm %s; rm %s; rm %s'%(outpath+'lepoint?.fits',outpath+'mepoint?.fits',outpath+'hepoint?.fits'))

hd=pf.open(outpath+'mepoint0_all.fits')
tm=hd[1].data.field(0)
tm_s,order = np.unique(tm,return_index=True)
yp=np.in1d(tm_s,atm_s)
hd[1].data=hd[1].data[order][yp]
hd[1].data=hd[1].data[hd[1].data.field(1)>-1000]
hd.writeto(outpath+'mepoint0.fits',overwrite=True)
hd.close();del hd


hd=pf.open(outpath+'mepoint1_all.fits')
tm=hd[1].data.field(0)
tm_s,order = np.unique(tm,return_index=True)
yp=np.in1d(tm_s,atm_s)
hd[1].data=hd[1].data[order][yp]
hd[1].data=hd[1].data[hd[1].data.field(1)>-1000]
hd.writeto(outpath+'mepoint1.fits',overwrite=True)
hd.close();del hd

hd=pf.open(outpath+'mepoint2_all.fits')
tm=hd[1].data.field(0)
tm_s,order = np.unique(tm,return_index=True)
yp=np.in1d(tm_s,atm_s)
hd[1].data=hd[1].data[order][yp]
hd[1].data=hd[1].data[hd[1].data.field(1)>-1000]
hd.writeto(outpath+'mepoint2.fits',overwrite=True)
hd.close();del hd

hdatt = pf.open(outpath+'pointatt.fits')
attm = hdatt[1].data.field(0)
yep = np.in1d(attm,tm_s)

hdatt[1].data = hdatt[1].data[yep]
hdatt[2].data = hdatt[2].data[yep]
hdatt[3].data = hdatt[3].data[yep]
hdatt[4].data = hdatt[4].data[yep]
hdatt.writeto(outpath+'pointatt.fits',overwrite=True)
hdatt.close();del hdatt


float('l')
hd0=pf.open(outpath+'mepoint0.fits')
hd1=pf.open(outpath+'mepoint1.fits')
hd2=pf.open(outpath+'mepoint2.fits')
tm0=hd0[1].data.field(0)
tm1=hd1[1].data.field(0)
tm2=hd2[1].data.field(0)

ct0=hd0[1].data.field(1)
ct1=hd1[1].data.field(1)
ct2=hd2[1].data.field(1)
er0=hd0[1].data.field(2)
er1=hd1[1].data.field(2)
er2=hd2[1].data.field(2)
cts = ct0+ct1+ct2




yp = (cts>3000)



rat = 3000./ct0[yp]
#norm = np.load('scox1box0.npy')
#rat = norm/ct0[yp]
tms = tm0[yp]

ayp = np.in1d(attm,tms)
hdatt = pf.open(outpath+'pointatt.fits')
ra = hdatt[1].data['Ra'][ayp]
dec= hdatt[1].data['Dec'][ayp]
hdatt.close();del hdatt

dui = [(ra[0],dec[0])]
duid = [[0]]
obra,obdec = 244.979,-15.640
for i in range(1,len(ra)):
    print i
    tpra = ra[i]
    tpdec = dec[i]
    flag=0
    for j in range(len(dui)):
        tpds = distance(tpra,tpdec,dui[j][0],dui[j][1])
        if tpds<0.002:
            duid[j].append(i)
            flag=1
            break
    if not flag:
        dui.append((tpra,tpdec))
        duid.append([i])

s=0;m=0
for ind,i in enumerate(duid):
    if len(i)>50:
        s+=1
    m = max(m,len(i))

print s,m
for ind,i in enumerate(duid):
    if len(i)==m:
        s=ind
        break

with open('/sharefs/hbkg/data/SCAN/nangyi/scox1_3/stepinfo.txt','w')as f:
    for i in duid[s]:
        print>>f,i


bct0 = ct0[yp]*rat
bct1 = ct1[yp]*rat
bct2 = ct2[yp]*rat
ber0 = er0[yp]*rat
ber1 = er1[yp]*rat
ber2 = er2[yp]*rat

col1 = pf.Column(name = 'Time', format='D',array = tms)
col2 = pf.Column(name = 'Counts', format='D',array = bct0)
col3 = pf.Column(name = 'Stat_err', format='D',array = ber0)
cols = pf.ColDefs([col1,col2,col3])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(outpath+'mepoint0u.fits',overwrite=True)

col1 = pf.Column(name = 'Time', format='D',array = tms)
col2 = pf.Column(name = 'Counts', format='D',array = bct1)
col3 = pf.Column(name = 'Stat_err', format='D',array = ber1)
cols = pf.ColDefs([col1,col2,col3])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(outpath+'mepoint1u.fits',overwrite=True)

col1 = pf.Column(name = 'Time', format='D',array = tms)
col2 = pf.Column(name = 'Counts', format='D',array = bct2)
col3 = pf.Column(name = 'Stat_err', format='D',array = ber2)
cols = pf.ColDefs([col1,col2,col3])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(outpath+'mepoint2u.fits',overwrite=True)



'''
hedatas0 =glob.glob(hepath+'P0101295*hebox0.fits')
hedatas0=sorted([i.split('/')[-1][:13] for i in hedatas0])

hedatas1 =glob.glob(hepath+'P0101295*hebox1.fits')
hedatas1=sorted([i.split('/')[-1][:13] for i in hedatas1])

hedatas2 =glob.glob(hepath+'P0101295*hebox2.fits')
hedatas2=sorted([i.split('/')[-1][:13] for i in hedatas2])

unidatas = sorted([i for i in set(hedatas0+hedatas1+hedatas2)])
outpath = '/sharefs/hbkg/data/SCAN/PSF/ME/295Low/'

#os.system('rm %s; rm %s; rm %s'%(outpath+'lepoint?_all.fits',outpath+'hepoint?_all.fits',outpath+'hepoint?_all.fits'))
hedatas = [hepath+j+'hebox0.fits' for j in hedatas0]
fitstogetherMulti(hedatas,outpath+'hepoint0_all.fits')

hedatas = [hepath+j+'hebox1.fits' for j in hedatas1]
fitstogetherMulti(hedatas,outpath+'hepoint1_all.fits')

hedatas = [hepath+j+'hebox2.fits' for j in hedatas2]
fitstogetherMulti(hedatas,outpath+'hepoint2_all.fits')


#atts = ['/sharefs/hbkg/data/SCAN/nangyi/ME_data/att/'+j+'crtatt' for j in unidatas]
#os.system('rm %s'%(outpath+'pointatt.fits'))
#fitstogetherMulti(atts,outpath+'pointatt.fits')
hdatt = pf.open(outpath+'pointatt.fits')
attm = hdatt[1].data.field(0)
hdatt.close();del hdatt

hd=pf.open(outpath+'hepoint0_all.fits')
tm=hd[1].data.field(0)
yp=np.in1d(tm,attm)
hd[1].data=hd[1].data[yp]
hd.writeto(outpath+'hepoint0.fits')
hd.close();del hd


hd=pf.open(outpath+'hepoint1_all.fits')
tm=hd[1].data.field(0)
yp=np.in1d(tm,attm)
hd[1].data=hd[1].data[yp]
hd.writeto(outpath+'hepoint1.fits')
hd.close();del hd

hd=pf.open(outpath+'hepoint2_all.fits')
tm=hd[1].data.field(0)
yp=np.in1d(tm,attm)
hd[1].data=hd[1].data[yp]
hd.writeto(outpath+'hepoint2.fits')
hd.close();del hd
'''
#os.system('rm %spointinfonew.txt; rm %s?epoint?.???'%(outpath,outpath))
