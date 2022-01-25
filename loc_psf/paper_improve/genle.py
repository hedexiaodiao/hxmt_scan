import numpy as np
import sys
import os
import glob
from astropy.io import fits as pf
from tool import *
def mkdir_try(dirname):
    if os.path.exists(dirname) == 0:
        try:
            os.makedirs(dirname)
        except OSError:
            print('Wrong with make dir:\n'+dirname)

if len(sys.argv)==3:
    energy_str = sys.argv[2]
    if energy_str=='0':
        print('use energy 1-6 keV')
        minpi =106 ### energy = 1 keV
        maxpi = 710  ### used by cwang
        outpath = '/sharefs/hbkg/user/luoqi/psfl/genlc/1_6/'
    elif energy_str=='1':
        print('use energy 2-6 keV')
        minpi = 225  ### energy = 2 keV
        maxpi = 710  ### used by cwang
        outpath = '/sharefs/hbkg/user/luoqi/psfl/genlc/2_6/'
    else:
        print('wrong input!')
else:
    print('use energy 1-6 keV')
    minpi = 106  ### energy = 1 keV
    maxpi = 710  ### used by cwang
    outpath = '/sharefs/hbkg/user/luoqi/psfl/genlc/1_6/'

obsid = sys.argv[1]
pfile = "/sharefs/hbkg/user/luoqi/HXMT_SCAN/loc_psf/pfiles/le%s"%str(obsid)
os.system('mkdir %s'%pfile);os.system("rm %s/*"%pfile);
pfilepath = "%s;/home/hxmt/hxmtsoft2/hxmtsoftv2.02/install/x86_64-pc-linux-gnu-libc2.12/syspfiles"%pfile
###cmdinit = 'sleep 1;source /sharefs/hbkg/user/nangyi/hxmtsoft_v2.02.sh ;export PFILES="%s"'%(pfilepath)
cmdinit = 'sleep 1;source /sharefs/hbkg/user/luoqi/home/hxmtsoft_v2.02.sh ;export PFILES="%s"'%(pfilepath)
###minpi = 226


basedatapath = '/hxmt/work/HXMT-DATA/1L/'+"A%s/%s/%s/"%(obsid[1:3],obsid[:8],obsid[:11])
datapath = glob.glob(basedatapath+obsid+'-*')[0]
basefile = get_rawdata(datapath, instrument="LE")

###outpath = '/sharefs/hbkg/user/luoqi/psfl/genlc/1_6/'
###outpath = '/sharefs/hbkg/user/luoqi/psfl/genlc/2_6/'
##outpath = '/sharefs/hbkg/data/SCAN/luoqi/LE_data/'
pipath = outpath + 'pi'
grdpath = outpath + 'grd'
gtipath = outpath + 'gti'
screenpath = outpath + 'screen'
lcpath = outpath + 'lc'
mkdir_try(pipath)
mkdir_try(grdpath)
mkdir_try(gtipath)
mkdir_try(screenpath)
mkdir_try(lcpath)

pifile = outpath + 'pi/%slepi.fits'%obsid
grdfile = outpath + 'grd/%slegrd.fits'%obsid
gtifile = outpath + 'gti/%slegti.fits'%obsid
newgtifile = outpath + 'gti/%slegtinew.fits'%obsid
screenfile = outpath + 'screen/%slescreen.fits'%obsid
lcfile = outpath + 'lc/%sle'%obsid
lcnametxt = outpath + 'lc/%slelc.txt'%obsid
bkgfile = outpath + 'lc/%slebkg'%obsid
#basefile['EHK']='/sharefs/hbkg/data/SCAN/nangyi/HE_point/ehk/%sehk.fits'%obsid

cmd=''
if not os.path.exists(pifile):
    cmd=' ; lepical evtfile=%s tempfile=%s outfile=%s'%(basefile['EVT'],basefile['TH'],pifile)
    #os.system(cmd)

if not os.path.exists(grdfile):
    cmd +=' ;lerecon evtfile=%s outfile=%s instatusfile=%s'%(pifile,grdfile,basefile['InsStat'])
    #os.system(cmd)
os.system(cmdinit + cmd)

start_time = 100000000
stop_time  = 400000000
gticmd=''
if 1:#not os.path.exists(gtifile):
    #cmd +=' ;legtigen evtfile=%s tempfile="%s" instatusfile="%s" ehkfile="%s" outfile="%s" defaultexpr=NONE expr="ELV>10&&DYE_ELV>10&&COR>8&&SAA_FLAG==0&&T_SAA>=200&&TN_SAA>=200&&ANG_DIST<=0.05" clobber=yes history=yes'%(basefile['EVT'],basefile['TH'],basefile['InsStat'], basefile['EHK'], gtifile)
    gticmd +=' ;legtigen evtfile=%s tempfile="%s" instatusfile="%s" ehkfile="%s" outfile="%s" defaultexpr=NONE expr="ELV>5" clobber=yes history=yes'%(basefile['EVT'],basefile['TH'],basefile['InsStat'], basefile['EHK'], gtifile)
    os.system(cmdinit + gticmd)
    hdul = pf.open(gtifile)
    data = hdul[1].data
    data['START'] = start_time
    data['STOP'] = stop_time
    hdul.writeto(gtifile, clobber=True)
    hdul.close()
    newgticmd = ' ;legti %s %s %s' % (grdfile, gtifile, newgtifile)
    os.system(cmdinit + newgticmd)

screencmd=''
if 1:#not os.path.exists(screenfile):
    screencmd +=' ;lescreen evtfile=%s gtifile=%s outfile=%s userdetid="0-95"'%(grdfile,gtifile,screenfile)
    os.system(cmdinit + screencmd)


thd = pf.open(screenfile)#,endcard=False)
try:
    print len(thd)
except:
    pass

startm = int(thd[2].data.field(0)[0])+1
stoptm = int(thd[2].data.field(1)[-1])-1

thd.close();del thd
###29,54,87###bad detector, 54 used in old
cmd = ';lelcgen evtfile=%s outfile=%s binsize=1 starttime=%s stoptime=%s minPI=%s maxPI=%s userdetid="0,2,3,4,6,7,8,9,10,12,14,20,22,23,24,25,26,28,30;32,34,35,36,38,39,40,41,42,44,46,52,55,56,57,58,60,61,62;64,66,67,68,70,71,72,73,74,76,78,84,86,88,89,90,92,93,94" eventtype=1'%(screenfile,lcfile,startm,stoptm,minpi,maxpi)
cmd += ' ; ls -t %s_g0*.lc | head -1 >  %s'%(lcfile,lcnametxt)
cmd += ' ; python lebkgmap_202005_box0.py lc %s %s %s %s %s %s_box0'%(screenfile,newgtifile,lcnametxt,minpi,maxpi,bkgfile)
cmd += ' ; python lebkgmap_202005_box1.py lc %s %s %s %s %s %s_box1'%(screenfile,newgtifile,lcnametxt,minpi,maxpi,bkgfile)
cmd += ' ; python lebkgmap_202005_box2.py lc %s %s %s %s %s %s_box2'%(screenfile,newgtifile,lcnametxt,minpi,maxpi,bkgfile)
os.system(cmdinit + cmd)
'''
ljygti = glob.glob('/sharefs/hbkg/user/liaojy/HXMT_BKG_LE/Data_filter_Auto/%s*/%s*_LE_GTI_Auto.txt'%(obsid,obsid))[0]
ljygti = np.loadtxt(ljygti)
ljytm = []
for i in ljygti:
    ljytm = np.append(ljytm,np.arange(i[0],i[1]))
'''
hd1=pf.open('%s_g0_0-30.lc'%lcfile)
hd2=pf.open('%s_box0.lc'%bkgfile)
orgtm = hd1[1].data['TIME']
yp = (hd1[1].data['FRACEXP']==1.)#&(np.in1d(orgtm,ljytm)))
tm = hd1[1].data['TIME'][yp]
ct = hd1[1].data['COUNTS'][yp]
cter = hd1[1].data['ERROR'][yp]
bkgct = hd2[1].data['RATE'][yp]
bkger = hd2[1].data['Error'][yp]
netct = ct-bkgct
neter = np.sqrt(cter**2+bkger**2)
col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = netct)
col3 = pf.Column(name = 'Stat_err', format='D',array = neter)
col4 = pf.Column(name = 'Cts_org', format='D',array = ct)
col5 = pf.Column(name = 'Cts_err', format='D',array = cter)
col6 = pf.Column(name = 'Bkg', format='D',array = bkgct)
col7 = pf.Column(name = 'Bkg_err', format='D',array = bkger)
cols = pf.ColDefs([col1,col2,col3,col4,col5,col6,col7])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box0.fits",overwrite=True)
hd1.close();hd2.close();del hd1,hd2


hd1=pf.open('%s_g1_32-62.lc'%lcfile)
hd2=pf.open('%s_box1.lc'%bkgfile)
orgtm = hd1[1].data['TIME']
yp = (hd1[1].data['FRACEXP']==1.)#&(np.in1d(orgtm,ljytm)))
tm = hd1[1].data['TIME'][yp]
ct = hd1[1].data['COUNTS'][yp]
cter = hd1[1].data['ERROR'][yp]
bkgct = hd2[1].data['RATE'][yp]
bkger = hd2[1].data['Error'][yp]
netct = ct-bkgct
neter = np.sqrt(cter**2+bkger**2)
col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = netct)
col3 = pf.Column(name = 'Stat_err', format='D',array = neter)
col4 = pf.Column(name = 'Cts_org', format='D',array = ct)
col5 = pf.Column(name = 'Cts_err', format='D',array = cter)
col6 = pf.Column(name = 'Bkg', format='D',array = bkgct)
col7 = pf.Column(name = 'Bkg_err', format='D',array = bkger)
cols = pf.ColDefs([col1,col2,col3,col4,col5,col6,col7])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box1.fits",overwrite=True)
hd1.close();hd2.close();del hd1,hd2

hd1=pf.open('%s_g2_64-94.lc'%lcfile)
hd2=pf.open('%s_box2.lc'%bkgfile)
orgtm = hd1[1].data['TIME']
yp = (hd1[1].data['FRACEXP']==1.)#&(np.in1d(orgtm,ljytm)))
tm = hd1[1].data['TIME'][yp]
ct = hd1[1].data['COUNTS'][yp]
cter = hd1[1].data['ERROR'][yp]
bkgct = hd2[1].data['RATE'][yp]
bkger = hd2[1].data['Error'][yp]
netct = ct-bkgct
neter = np.sqrt(cter**2+bkger**2)
col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = netct)
col3 = pf.Column(name = 'Stat_err', format='D',array = neter)
col4 = pf.Column(name = 'Cts_org', format='D',array = ct)
col5 = pf.Column(name = 'Cts_err', format='D',array = cter)
col6 = pf.Column(name = 'Bkg', format='D',array = bkgct)
col7 = pf.Column(name = 'Bkg_err', format='D',array = bkger)
cols = pf.ColDefs([col1,col2,col3,col4,col5,col6,col7])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box2.fits",overwrite=True)
hd1.close();hd2.close();del hd1,hd2



