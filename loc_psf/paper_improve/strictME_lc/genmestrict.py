import sys
import numpy as np
import xlrd
import os
import glob
from astropy.io import fits as pf
#bdfile = '/home/hxmt/hxmtsoft2/hxmtsoftv2.02/install/x86_64-pc-linux-gnu-libc2.12/refdata/medetectorstatus.fits'
from tool import *
sys.path.append("/sharefs/hbkg/user/nangyi/HXMT_GPS_Tools")
from Lctools_202004 import correct
def mkdir_try(dirname):
    if os.path.exists(dirname) == 0:
        try:
            os.makedirs(dirname)
        except OSError:
            print('Wrong with make dir:\n'+dirname)

if len(sys.argv)>=3:
    energy_str = sys.argv[2]
    cut_value = sys.argv[3]
    if energy_str=='7_12':
        print('use energy 7-12 keV')
        minpi = 68  ###energy =7 keV
        maxpi = 154 ###energy = 12 keV
        outpath = '/sharefs/hbkg/user/luoqi/psfl/genlc_strictME/7_12_%s/'%str(cut_value)
    elif energy_str=='7_20':
        print('use energy 7-20 keV')
        minpi = 68  ###energy =7 keV
        maxpi = 290  ###energy = 20 keV
        outpath = '/sharefs/hbkg/user/luoqi/psfl/genlc_strictME/7_20_%s/'%str(cut_value)
    elif energy_str=='7_40':
        print('use energy 7-40 keV')
        minpi = 68
        maxpi = 631
        outpath = '/sharefs/hbkg/user/luoqi/psfl/genlc_strictME/7_40_%s/'%str(cut_value)
    else:
        print('wrong input!')
else:
    print('use energy 7-12 keV')
    minpi = 68  ###energy =7 keV
    maxpi = 154  ###energy = 12 keV
    outpath = '/sharefs/hbkg/user/luoqi/psfl/genlc/7_12/'

obsid = sys.argv[1]

pfile = "/sharefs/hbkg/user/luoqi/HXMT_SCAN/loc_psf/pfiles/me%s"%str(obsid)
os.system('mkdir %s'%pfile);os.system("rm %s/*"%pfile);
pfilepath = "%s;/home/hxmt/hxmtsoft2/hxmtsoftv2.02/install/x86_64-pc-linux-gnu-libc2.12/syspfiles"%pfile
###cmdinit = 'sleep 1;source /sharefs/hbkg/user/nangyi/hxmtsoft_v2.02.sh ;export PFILES="%s" '%(pfilepath)
cmdinit = 'sleep 1;source /sharefs/hbkg/user/luoqi/home/hxmtsoft_v2.02.sh ;export PFILES="%s"'%(pfilepath)

basedatapath = '/hxmt/work/HXMT-DATA/1L/'+"A%s/%s/%s/"%(obsid[1:3],obsid[:8],obsid[:11])
datapath = glob.glob(basedatapath+obsid+'-*')[0]
basefile = get_rawdata(datapath, instrument="ME")

###outpath = '/sharefs/hbkg/user/luoqi/psfl/genlc/7_12/'#'/sharefs/hbkg/data/SCAN/luoqi/ME_data/7_40/'
###outpath = '/sharefs/hbkg/user/luoqi/psfl/genlc/7_20/'
pipath = outpath + 'pi'
grdpath = outpath + 'grd'
gtipath = outpath + 'gti'
screenpath = outpath + 'screen'
lcpath = outpath + 'lc'
attpath = outpath +'att'
specpath = outpath + 'spec'
mkdir_try(pipath)
mkdir_try(grdpath)
mkdir_try(gtipath)
mkdir_try(screenpath)
mkdir_try(lcpath)
mkdir_try(attpath)
mkdir_try(specpath)

pifile = outpath + 'pi/%smepi.fits'%obsid
grdfile = outpath + 'grd/%smegrd.fits'%obsid
deadfile = outpath + 'grd/%smedead.fits'%obsid
gtifile = outpath + 'gti/%smegti.fits'%obsid
newgtifile = outpath + 'gti/%smegtinew.fits'%obsid
newbdfile = outpath + 'gti/%smebd.fits'%obsid
screenfile = outpath + 'screen/%smescreen.fits'%obsid
specfile = outpath + 'spec/%smespec'%obsid
specnametxt =  outpath + 'spec/%smespec.txt'%obsid
specbkgfile = outpath + 'spec/%sbkgspec'%obsid
rspfile = outpath + 'spec/%srsp'%obsid
lcfile = outpath + 'lc/%sme'%obsid
lcnametxt = outpath + 'lc/%smelc.txt'%obsid
bkgfile = outpath + 'lc/%smebkg'%obsid
#basefile['EHK']='/sharefs/hbkg/data/SCAN/nangyi/HE_point/ehk/%sehk.fits'%obsid

crtattfile = outpath + 'att/%scrtatt'%obsid
correct(basefile['ATT'],crtattfile)

cmd=''
if not os.path.exists(pifile):
    cmd='; mepical evtfile=%s tempfile=%s outfile=%s'%(basefile['EVT'],basefile['TH'],pifile)

if not os.path.exists(grdfile):
    cmd +=' ;megrade evtfile=%s outfile=%s deadfile=%s binsize=1'%(pifile,grdfile,deadfile)
os.system(cmdinit + cmd)

start_time = 100000000
stop_time  = 400000000
gticmd=''
if 1:#not os.path.exists(gtifile):
    #cmd +=' ;megtigen tempfile="%s" ehkfile="%s" outfile="%s" defaultexpr=NONE expr="ELV>8&&COR>8&&SAA_FLAG==0&&T_SAA>300&&TN_SAA>300&&ANG_DIST<=0.02" clobber=yes history=yes'%(basefile['TH'], basefile['EHK'], gtifile)
    gticmd +=' ;megtigen tempfile="%s" ehkfile="%s" outfile="%s" defaultexpr=NONE expr="ELV>8" clobber=yes history=yes'%(basefile['TH'], basefile['EHK'], gtifile)
    os.system(cmdinit + gticmd)
    hdul = pf.open(gtifile)
    data = hdul[1].data
    data['START'] = start_time
    data['STOP'] = stop_time
    hdul.writeto(gtifile, clobber=True)
    hdul.close()
    #newgticmd =' ;megti %s %s %s $HEADAS/refdata/medetectorstatus.fits %s'%(grdfile, gtifile, newgtifile, newbdfile)
    new_strictbdfile = outpath + 'gti/%s_ME_status_new.fits'%obsid
    newgticmd = ' ;python megti_detid.py %s %s %s $HEADAS/refdata/medetectorstatus.fits %s %s'%(grdfile, gtifile, newgtifile, new_strictbdfile, cut_value)
    os.system(cmdinit + newgticmd)

newbdfile = outpath + 'gti/mebadall.fits'

#new_strictbdfile = outpath + 'gti/%s_ME_status_new.fits'%obsid

screencmd=''
if not os.path.exists(screenfile):
    #screencmd +=' ;mescreen evtfile=%s gtifile=%s outfile=%s baddetfile=%s userdetid="0-53"'%(grdfile,gtifile,screenfile,newbdfile)
    screencmd +=' ;mescreen evtfile=%s gtifile=%s outfile=%s baddetfile=%s userdetid="0-53"'%(grdfile,gtifile,screenfile,new_strictbdfile)
print(screencmd)
os.system(cmdinit + screencmd)
#plt.switch_backend('tkagg')
'''
#newbdfile = outpath + "gti/mebadall.fits"
bdnum_r = genmebd(screenfile,outpath + 'gti/%smebdnew.fits'%obsid,newbdfile)
newbdfile = outpath + 'gti/%smebdnew.fits'%obsid
sys.exit(0)
#hd=pf.open(newbdfile)
#bdnum=hd[1].data.field(0)
#hd[1].data = hd[1].data[bdnum>=576]
#hd.writeto(newbdfile,overwrite=True)
#hd.close();del hd

#if not os.path.exists(screenfile):
cmd =' ;mescreen evtfile=%s gtifile=%s outfile=%s baddetfile=%s userdetid="0-53"'%(grdfile,newgtifile,screenfile,newbdfile)

os.system(cmdinit + cmd)
'''
thd = pf.open(screenfile)#,endcard=False)
try:
    print len(thd)
except:
    pass

'''
start_list = [259271900,259271965,259271865,259271995]
stop_list = [259271940,259271995,259271900,259272020]
startm = start_list[3]
stoptm = stop_list[3]
'''
startm = int(thd[2].data.field(0)[0])+1
stoptm = int(thd[2].data.field(1)[-1])-1
#'''

thd.close();del thd
#float('l')
cmd = '; melcgen evtfile=%s deadfile=%s outfile=%s userdetid="0-7,11-17;18-25,29-35;36-43,47-53" starttime=%s stoptime=%s minPI=%s maxPI=%s binsize=1 deadcorr=yes'%(screenfile,deadfile,lcfile,startm,stoptm,minpi,maxpi)
cmd += ' ; ls -t %s_g0*.lc | head -1 >  %s'%(lcfile,lcnametxt)
cmd += ' ; python mebkgmap_2012005_lc_err.py lc %s %s %s %s %s %s %s %s %s %s'%(screenfile,basefile['EHK'],newgtifile,deadfile,basefile['TH'],lcnametxt,minpi,maxpi,bkgfile,new_strictbdfile)
print(cmd)
os.system(cmdinit + cmd)

burst_ra = 263.353
burst_dec = -33.389
speccmd = '; mespecgen evtfile=%s deadfile=%s outfile=%s userdetid="0-7,11-17;18-25,29-35;36-43,47-53" starttime=%s stoptime=%s minPI=%s maxPI=%s'%(screenfile,deadfile,specfile,startm,stoptm,minpi,maxpi)
speccmd += ' ; ls -t %s_g0*.pha | head -1 >  %s'%(specfile,specnametxt)
speccmd += ' ; python mebkgmap_2012005_lc_err.py spec %s %s %s %s %s %s %s %s %s %s'%(screenfile,basefile['EHK'],newgtifile,deadfile,basefile['TH'],specnametxt,minpi,maxpi,specbkgfile,new_strictbdfile)
speccmd += ' ; merspgen phafile=%s_g0_0-17.pha outfile=%s attfile=%s ra=%s dec=%s'%(specfile,rspfile,crtattfile,burst_ra,burst_dec)
print(speccmd)
os.system(cmdinit + speccmd)

###sfile=pf.open("%s"%newbdfile)
sfile=pf.open("%s"%new_strictbdfile)
bdnum=sfile[1].data.field(0)
bdbox0 = np.in1d(range(256)+range(352,576),bdnum).sum()
bdbox1 = np.in1d(range(576,832)+range(928,1152),bdnum).sum()
bdbox2 = np.in1d(range(1152,1408)+range(1504,1728),bdnum).sum()
sfile.close();del sfile
r0=480./(480-bdbox0)
r1=480./(480-bdbox1)
r2=480./(480-bdbox2)
rall=[r0,r1,r2]
print np.array(range(256)+range(352,576))[~np.in1d(range(256)+range(352,576),bdnum)],np.array(range(256)+range(352,576))[~np.in1d(range(256)+range(352,576),bdnum)].shape
###r0,r1,r2=1,1,1####don't know why use 1,1,1

hd1=pf.open('%s_g0_0-17.lc'%lcfile)
hd2=pf.open('%s_box0.lc'%bkgfile)
yp = (hd1[1].data['FRACEXP']==1.)
tm = hd1[1].data['TIME'][yp]
ct = hd1[1].data['COUNTS'][yp]*r0
cter = hd1[1].data['ERROR'][yp]*r0
bkgct = hd2[1].data['RATE'][yp]
bkger = hd2[1].data['Error'][yp]
print ct.sum()/bkgct.sum()
with open('bkg0.txt','a')as f:
    print>>f,ct.sum(),bkgct.sum(),obsid

netct = ct-bkgct
neter = np.sqrt(cter**2+bkger**2)
col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = netct)
col3 = pf.Column(name = 'Stat_err', format='D',array = neter)
#col4 = pf.Column(name = 'Cts_org', format='D',array = ct)
#col5 = pf.Column(name = 'Cts_err', format='D',array = cter)
#col6 = pf.Column(name = 'Bkg', format='D',array = bkgct)
#col7 = pf.Column(name = 'Bkg_err', format='D',array = bkger)
col4 = pf.Column(name = 'OGCounts', format='D',array = ct)
col5 = pf.Column(name = 'BGCounts', format='D',array = bkgct)

cols = pf.ColDefs([col1,col2,col3,col4,col5])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box0.fits",overwrite=True)
hd1.close();hd2.close();del hd1,hd2


hd1=pf.open('%s_g1_18-35.lc'%lcfile)
hd2=pf.open('%s_box1.lc'%bkgfile)
yp = (hd1[1].data['FRACEXP']==1.)
tm = hd1[1].data['TIME'][yp]
ct = hd1[1].data['COUNTS'][yp]*r1
cter = hd1[1].data['ERROR'][yp]*r1
bkgct = hd2[1].data['RATE'][yp]
bkger = hd2[1].data['Error'][yp]
netct = ct-bkgct
print ct.sum()/bkgct.sum()
with open('bkg1.txt','a')as f:
    print>>f,ct.sum(),bkgct.sum(),obsid


neter = np.sqrt(cter**2+bkger**2)
col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = netct)
col3 = pf.Column(name = 'Stat_err', format='D',array = neter)
#col4 = pf.Column(name = 'Cts_org', format='D',array = ct)
#col5 = pf.Column(name = 'Cts_err', format='D',array = cter)
#col6 = pf.Column(name = 'Bkg', format='D',array = bkgct)
#col7 = pf.Column(name = 'Bkg_err', format='D',array = bkger)
col4 = pf.Column(name = 'OGCounts', format='D',array = ct)
col5 = pf.Column(name = 'BGCounts', format='D',array = bkgct)

cols = pf.ColDefs([col1,col2,col3,col4,col5])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box1.fits",overwrite=True)
hd1.close();hd2.close();del hd1,hd2


hd1=pf.open('%s_g2_36-53.lc'%lcfile)
hd2=pf.open('%s_box2.lc'%bkgfile)
yp = (hd1[1].data['FRACEXP']==1.)
tm = hd1[1].data['TIME'][yp]
ct = hd1[1].data['COUNTS'][yp]*r2
cter = hd1[1].data['ERROR'][yp]*r2
bkgct = hd2[1].data['RATE'][yp]
bkger = hd2[1].data['Error'][yp]
netct = ct-bkgct
print ct.sum()/bkgct.sum()
with open('bkg2.txt','a')as f:
    print>>f,ct.sum(),bkgct.sum(),obsid


neter = np.sqrt(cter**2+bkger**2)
col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = netct)
col3 = pf.Column(name = 'Stat_err', format='D',array = neter)
#col4 = pf.Column(name = 'Cts_org', format='D',array = ct)
#col5 = pf.Column(name = 'Cts_err', format='D',array = cter)
#col6 = pf.Column(name = 'Bkg', format='D',array = bkgct)
#col7 = pf.Column(name = 'Bkg_err', format='D',array = bkger)
col4 = pf.Column(name = 'OGCounts', format='D',array = ct)
col5 = pf.Column(name = 'BGCounts', format='D',array = bkgct)

cols = pf.ColDefs([col1,col2,col3,col4,col5])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box2.fits",overwrite=True)
hd1.close();hd2.close();del hd1,hd2



