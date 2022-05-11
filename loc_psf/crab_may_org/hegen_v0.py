import numpy as np
import sys
import os
import glob
from astropy.io import fits as pf
from tool import *
sys.path.append("/sharefs/hbkg/user/nangyi/HXMT_GPS_Tools")
from Lctools_202004 import correct
def mkdir_try(dirname):
    if os.path.exists(dirname) == 0:
        try:
            os.makedirs(dirname)
        except OSError:
            print('Wrong with make dir:\n'+dirname)

obsid = sys.argv[1]
pfile = "/sharefs/hbkg/user/luoqi/HXMT_SCAN/loc_psf/pfilesv2/he%s"%str(obsid)
os.system('mkdir %s'%pfile);os.system("rm %s/*"%pfile);
pfilepath = "%s;/home/hxmt/hxmtsoft2/hxmtsoftv2.02/install/x86_64-pc-linux-gnu-libc2.12/syspfiles"%pfile
cmdinit = 'sleep 1;source /sharefs/hbkg/user/nangyi/hxmtsoft_v2.02.sh ;export PFILES="%s" '%(pfilepath)

basedatapath = '/hxmt/work/HXMT-DATA/1L/'+"A%s/%s/%s/"%(obsid[1:3],obsid[:8],obsid[:11])
datapath = glob.glob(basedatapath+obsid+'-*')[0]
basefile = get_rawdata(datapath, instrument="HE")
outpath = '/sharefs/hbkg/data/SCAN/luoqi/Crabv0_HE/'
pipath = outpath + 'pi'
grdpath = outpath + 'grd'
gtipath = outpath + 'gti'
screenpath = outpath + 'screen'
lcpath = outpath + 'lc'
attpath = outpath + 'att'
ehkpath = outpath + 'ehk'
mkdir_try(pipath)
mkdir_try(grdpath)
mkdir_try(gtipath)
mkdir_try(screenpath)
mkdir_try(lcpath)
mkdir_try(attpath)
mkdir_try(ehkpath)

pifile = outpath + 'pi/%shepi.fits'%obsid
gtifile = outpath + 'gti/%shegti.fits'%obsid
newgtifile = outpath + 'gti/%shegtinew.fits'%obsid
screenfile = outpath + 'screen/%shescreen.fits'%obsid
lcfile = outpath + 'lc/%she'%obsid
lcnametxt = outpath + 'lc/%shelc.txt'%obsid
bkgfile = outpath + 'lc/%shebkg'%obsid
crtattfile = outpath + 'att/%scrtatt'%obsid
if not os.path.exists(crtattfile):
    correct(basefile['ATT'],crtattfile)
cmd=''
if not os.path.exists(pifile):
    cmd='; hepical evtfile="%s" outfile="%s" clobber=yes'%(basefile['EVT'],pifile)

if not os.path.exists(gtifile):
    cmd += ' ; hegtigen hvfile="%s" tempfile="%s" pmfile="%s" outfile="%s" ehkfile="%s" defaultexpr=NONE expr="ELV>10&&COR>8&&SAA_FLAG==0&&TN_SAA>300&&T_SAA>300&&ANG_DIST<=0.05" pmexpr="" clobber=yes history=yes'%(basefile['HV'], basefile['TH'], basefile['PM'],gtifile, basefile['EHK'])

if not os.path.exists(screenfile):
    cmd += ' ; hescreen evtfile="%s" gtifile="%s" outfile="%s" userdetid="0-17" eventtype=1 anticoincidence=yes starttime=0 stoptime=0 minPI=0 maxPI=255 clobber=yes history=yes'%(pifile,gtifile,screenfile)

os.system(cmdinit + cmd)
thd = pf.open(screenfile)#,endcard=False)
try:
    print len(thd)
except:
    pass

startm = int(thd[2].data.field(0)[0])+1
stoptm = int(thd[2].data.field(1)[-1])-1

thd.close();del thd

#minpi = [22,23,22,23,22,22,23,22,22,22,22,22,22,22,22,22,23,22]
#maxpi =[80,85,81,85,83,83,84,82,82,83,82,82,82,82,80,82,83,82]
###minpi = [22, 23, 22, 23, 22, 22, 23, 22, 22, 22, 22, 22, 22, 22, 22, 22, 23, 22]
###maxpi = [80, 85, 81, 85, 83, 84, 87, 82, 82, 83, 82, 82, 82, 81, 80, 82, 83, 82]
minpi = [22,23,22,23,22,22,23,22,22,22,22,22,22,22,22,22,23,22]
maxpi = [80,85,81,85,83,83,84,82,82,83,82,82,82,82,80,82,83,82]
###minpi = [23]*18
###maxpi =[82]*18
cout,bkg = [],[]
cmd=';pwd'
for i in range(18):
    cmd += ' ; helcgen evtfile="%s" outfile="%s" deadfile="%s" userdetid="%s" binsize=1 starttime=%s stoptime=%s minPI=%s maxPI=%s deadcorr=yes clobber=yes'%(screenfile,lcfile,basefile['DTime'],i,startm,stoptm,minpi[i],maxpi[i])

os.system(cmdinit + cmd)

for i in range(18):
    hd1=pf.open('%s_g0_%s.lc'%(lcfile,i))
    yp = (hd1[1].data['FRACEXP']==1)
    tm = hd1[1].data['TIME'][yp]
    ct = hd1[1].data['COUNTS'][yp]
    cter = hd1[1].data['ERROR'][yp]
    cout.append([i,tm,ct,cter])
    hd1.close()
    del hd1

cmd=';ls -t %s_g0_0.lc | head -1 >  %s'%(lcfile,lcnametxt)
cmd+=' ; python hebkgmap_id14_202005.py lc %s %s %s %s %s %s %s %s'%(screenfile, basefile['EHK'], gtifile,basefile['DTime'],lcnametxt,minpi[i], maxpi[i], bkgfile)
os.system(cmdinit + cmd)
for i in range(18):
    if i==16:continue;
    hd2=pf.open('%s_ID%s.lc'%(bkgfile,i))
    bgct = hd2[1].data['Rate'][yp]
    bgcter = hd2[1].data['Error'][yp]
    bkg.append([i,tm,bgct,bgcter])
    hd2.close()
    del hd2

cout=np.array(cout);bkg=np.array(bkg)
ct1 = cout[np.in1d(cout[:,0],[0,4,7,14,17])][:,2].sum()
er1 = np.sqrt((cout[np.in1d(cout[:,0],[0,4,7,14,17])][:,3]**2).sum())
ct2 = cout[np.in1d(cout[:,0],[1,3,8,11,12])][:,2].sum()
er2 = np.sqrt((cout[np.in1d(cout[:,0],[1,3,8,11,12])][:,3]**2).sum())
ct3 = cout[np.in1d(cout[:,0],[5,6,10,13,15])][:,2].sum()
er3 = np.sqrt((cout[np.in1d(cout[:,0],[5,6,10,13,15])][:,3]**2).sum())

bgct1 = bkg[np.in1d(bkg[:,0],[0,4,7,14,17])][:,2].sum()
bger1 = np.sqrt((bkg[np.in1d(bkg[:,0],[0,4,7,14,17])][:,3]**2).sum())
bgct2 = bkg[np.in1d(bkg[:,0],[1,3,8,11,12])][:,2].sum()
bger2 = np.sqrt((bkg[np.in1d(bkg[:,0],[1,3,8,11,12])][:,3]**2).sum())
bgct3 = bkg[np.in1d(bkg[:,0],[5,6,10,13,15])][:,2].sum()
bger3 = np.sqrt((bkg[np.in1d(bkg[:,0],[5,6,10,13,15])][:,3]**2).sum())

col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = ct1-bgct1)
col3 = pf.Column(name = 'Stat_err', format='D',array = np.sqrt(er1**2+bger1**2))
col4 = pf.Column(name = 'Cts_org', format='D',array = ct1)
col5 = pf.Column(name = 'Cts_err', format='D',array = er1)
col6 = pf.Column(name = 'Bkg', format='D',array = bgct1)
col7 = pf.Column(name = 'Bkg_err', format='D',array = bger1)
cols = pf.ColDefs([col1,col2,col3,col4,col5,col6,col7])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box0.fits",overwrite=True)

col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = ct2-bgct2)
col3 = pf.Column(name = 'Stat_err', format='D',array = np.sqrt(er2**2+bger2**2))
col4 = pf.Column(name = 'Cts_org', format='D',array = ct2)
col5 = pf.Column(name = 'Cts_err', format='D',array = er2)
col6 = pf.Column(name = 'Bkg', format='D',array = bgct2)
col7 = pf.Column(name = 'Bkg_err', format='D',array = bger2)
cols = pf.ColDefs([col1,col2,col3,col4,col5,col6,col7])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box1.fits",overwrite=True)

col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = ct3-bgct3)
col3 = pf.Column(name = 'Stat_err', format='D',array = np.sqrt(er3**2+bger3**2))
col4 = pf.Column(name = 'Cts_org', format='D',array = ct3)
col5 = pf.Column(name = 'Cts_err', format='D',array = er3)
col6 = pf.Column(name = 'Bkg', format='D',array = bgct3)
col7 = pf.Column(name = 'Bkg_err', format='D',array = bger3)
cols = pf.ColDefs([col1,col2,col3,col4,col5,col6,col7])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box2.fits",overwrite=True)

