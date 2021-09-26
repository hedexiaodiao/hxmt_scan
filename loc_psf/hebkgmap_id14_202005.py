#!/usr/bin/env python
#!coding=utf-8
'''
Model constructed by Background Group.
Lian Jinyuan, Zhangshu, Guo Chengcheng, Jin Jing, Zhangjuan, Zhang Shu, et al.
Mail liaojinyuan@ihep.ac.cn

'''
'''
This version was written by Ge Mingyu
Mail gemy@ihep.ac.cn

Usage:

hebkgmap lc/spec blind_det.FITS ehkfile.fits gtifile.fits deadtime.fits lcname/specname outnam_prefix (spec_time_arnge)
    lc/spec: lc for background lightcurve and spec for background light curve 
    screen.FITS: should only include the events for blind detecters.
    ehkfile.fits: the EHK file for the observation
    gtifile.fits: the GTI file for ME
    deadtime.fits: the Dead Time for ME
    lcname/specname: is ASCII file, which includes the name of the source file for small FOV
    chmin: minimum channel for light curve
    chmax: maximum channel for light curve
    outnam_prefix: the output prefix for the spectrum

Using interactive method in prompt.	

2019-04-10:
Add new paramter input
Upadte the lightcurve generation.
Upadte the error estimation.


2019-05-07:
Upadte data base for the model: 
(1) map data for 72*18
(2) decay components for  longer term
(3) Channel relation for 17 detectors with the blind detector
Update the procedure for background estimation:
(1) estimate the spectrum of the blind detector according the backgroud model (basic spectrum + the longer term variation) for specific the grid
(2) the real spetrum for this grid
(3) calculate the ratio between the estimated spectrum and real spectrum for 256 channels
(4) obtain the model spectrum for specific detector similar with the blind detector: spec_mod
(5) correct the 'spec_mod' according to the energy relation between one detector and the blind detector

2019-05-09:
Channel 10-20: set 0

2019-05-30:
Errors:
The error for spectrum is wrong. Form spec*r --> specerr*r[sindex] + spec*rerr[sindex]


2019-07-15:
Upate: the background for PWD: 54--70

2019-10-01:
Upate: Add corrections for lightcurves
2019-11-11:
Add lightcurve smoothing
Correct the error for light_curce calculation '_cor'

2019-12-10:
Remove the cnts==0 for Blind detector

2019-12-12:
overwrite deleted

2020-01-13:
Test for the GitLab

2020-04-29:
Change rr0 calculation:
        cnt0 = (np.sum(tmpspec0[220:256]))
        cnt1 = (np.sum(tmpspec1[220:256]))
        if iidd==14:
            cnt0 = (np.sum(tmpspec0[220:253]))
            cnt1 = (np.sum(tmpspec1[220:253]))

2020-05-13:

Add the background lightcurves of all detectors

'''

from astropy.io import fits as pf
import numpy as np
import matplotlib.pyplot as plt
import os
#import commands
import sys
import time
from scipy import interpolate
from xml.etree.ElementTree import ElementTree,Element

try:
    # Python 2
    xrange
except NameError:
    # Python 3, xrange is now named range
    xrange = range
    raw_input = input

SVer = '2.02'
BVer = '2.0.10'

print( "*********************************************************" )
print( "******************  Running HXMT Bkg   ******************" )
print( "*********************************************************" )
print( "*********************************************************" )
print( "*********************************************************" )
print( "************ PRINT: hebkgmap -h for usage   *************" )
print( "*********************************************************" )
print( "HXMTDAS version ", SVer," HXMT background for Insight-HXMT/HE, ver-", BVer )
print( "The energy range for background lightcurve should be the same as source lightcurve")

uage_method1 = 'Method 1: hebkgmap lc/spec screen.FITS ehkfile.fits gtifile.fits deadtime.fits lcname/specname chmin chmax outnam_prefix'
uage_method2 = 'Method 2: Using interactive method in prompt.'
uage_method3 = 'Method 3: hebkgmap sflag=lc/spec evtfile=screen.FITS ehkfile=ehkfile.fits gtifile=gtifile.fits dtname=deadtime.fits srcdat=lcname/specname chmin=chmin chmax=chmax outnam=outnam_prefix'

def print_usage(uage_method1,uage_method2,uage_method3):
    print(uage_method1)
    print(uage_method2)
    print(uage_method3)

def check_argument():
    sp_lc_select  = []
    evtfilename   = []
    ehkname       = []
    gtifile       = []
    dtname        = []
    sl_name       = []
    chmin         = []
    chmax         = []
    outnam        = []
    slgti         = []
    len_arg = len(sys.argv)
    if len_arg <=1:
        raise IOError("Error input argument, RUN 'hebkgmap -h' for help")
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print_usage(uage_method1,uage_method2)
        return False, False, False, False
    if((len_arg>1)&(len_arg<9)):
        raise IOError("Error input argument, RUN 'hebkgmap -h' for help")
        sys.exit()
    sysflag = 0
    for i in xrange(len_arg):
        arg = sys.argv[i]
        if i == 0:continue
        arg_split = arg.split('=')
        if (len(arg_split)==1):continue
        argname = arg_split[0]
        argval  = arg_split[1]
        if argname == 'sflag': 
            sp_lc_select.append(argval)
            sysflag = sysflag + 1
        if argname == 'evtfile': 
            evtfilename.append(argval)
            sysflag = sysflag + 1
        if argname == 'ehkfile':
            ehkname.append( argval )
            sysflag = sysflag + 1
        if argname == 'gtifile':
            gtifile.append(argval)
            sysflag = sysflag + 1
        if argname == 'deadtime':
            dtname.append(argval)
            sysflag = sysflag + 1
        if argname == 'srcdat':
            sl_name.append(argval.strip())
            sysflag = sysflag + 1
        if argname == 'chmin':
            chmin.append(int(argval))
            sysflag = sysflag + 1
        if argname == 'chmax':
            chmax.append(int(argval))
            sysflag = sysflag + 1
        if argname == 'outnam':
            outnam.append(argval)
            sysflag = sysflag + 1
        if (len(sys.argv)==10):
            slgti.append('')
        '''if (len(sys.argv)==11):
            if argname == 'newgti':
                slgti.append(argval)
                sysflag = sysflag + 1
        '''
    if((sysflag>0)&(sysflag<9)):
        raise IOError("Error input argument, RUN 'hebkgmap -h' for help")
        sys.exit()
    if ((sysflag==0) & (len_arg>=10)):
        sp_lc_select.append(sys.argv[1])
        evtfilename.append(sys.argv[2])
        ehkname.append(sys.argv[3])
        gtifile.append(sys.argv[4])
        dtname.append(sys.argv[5])
        sl_name.append(sys.argv[6])
        chmin.append(int(sys.argv[7]))
        chmax.append(int(sys.argv[8]))
        outnam.append(sys.argv[9])
        slgti.append('')
        '''if (len(sys.argv)==10):
            slgti.append('')
        if (len(sys.argv)==11):
            slgti.append(sys.argv[10])
        '''
    return sp_lc_select[0],evtfilename[0],ehkname[0],gtifile[0],dtname[0],sl_name[0],chmin[0],chmax[0],outnam[0],slgti[0]



if len(sys.argv)==2:
    if sys.argv[1]=='-h':
        print_usage(uage_method1,uage_method2,uage_method3)
    sys.exit()
elif len(sys.argv)>=2:
    sp_lc_select,evtfilename,ehkname,gtifile,dtname,sl_name,chmin,chmax,outnam,slgti=check_argument()
else:
    sp_lc_select= str(raw_input("Selection(spec/lc):"))
    evtfilename = str(raw_input("Screened events file:"))
    ehkname     = str(raw_input("EHK file:"))
    gtifile     = str(raw_input("GTI file:"))
    dtname      = str(raw_input("Dead time correction file:"))
    sl_name     = str(raw_input("FileName inluding spectra or lightcurve name:"))
    chmin       = int(raw_input("Minimum channel:"))
    chmax       = int(raw_input("Maximum channel:"))
    outnam      = str(raw_input("The prefix of output file name:"))
    #slgti       = str(raw_input("Specific time range file(NONE):"))
    slgti.append('')

HEADAS=os.getenv('HEADAS')
#HEADAS=os.getenv('REFPATH')
REFPATH=HEADAS+'/refdata/'

xml = ElementTree()
xml.parse(sys.argv[10])
chmin = [int(i) for i in xml.find('PARAMETER/MinPiPara').text.split(',')]
chmax = [int(i) for i in xml.find('PARAMETER/MaxPiPara').text.split(',')]

#chmin,chmax = [26, 27, 26, 27, 26, 27, 28, 26, 26, 26, 26, 26, 26, 26, 26, 26, 27, 26],[80,85,81,85,83,83,84,82,82,83,82,82,82,82,80,82,83,82]

print '*********************************************'
print chmin
print chmax
print '*********************************************'

'''Check the time range in GTI'''
def is_ingti(START,STOP,tl,tu):
    num = len(START)
    is_gti = 0
    for ii in xrange(0,num):
        t0=START[ii]
        t1=STOP[ii]
        flag0 = (tl<=t0) & (tu>=t0)
        flag1 = (tl>=t0) & (tu<=t1)
        flag2 = (tl<=t1) & (tu>=t1)
        if flag0 | flag1 | flag2:
            is_gti = 1
            return is_gti
    return is_gti
'''Check the time in GTI '''
def is_ingti2(START,STOP,ti):
    num = np.size(START)
    is_gti = 0
    if num == 1:
        t0=START
        t1=STOP
        flag0 = (ti>=t0) & (ti<=t1)
        if flag0:
            is_gti = 1
            return is_gti
    if num >= 2:
        for ii in xrange(0,num):
            t0=START[ii]
            t1=STOP[ii]
            flag0 = (ti>=t0) & (ti<=t1)
            if flag0:
                is_gti = 1
                return is_gti
    return is_gti
'''Give the tag for specific time'''
def time_gtiflag(time,START,STOP,tflag):
    for ii in xrange(0,len(time)):
        tflag[ii] = is_ingti2(START,STOP,time[ii])

'''Give the index for specific time'''
def flag_selection(tflag,tindex):
    cnt = 0
    for ii in xrange(0,len(tflag)):
        if (tflag[ii] == 1):
            tindex[cnt] = ii
            cnt = cnt+1
    return cnt

'''Give the index for bkg map'''
def bkgmap_index(LON,LAT,ascend_flag ,step):
    lon_index = int(LON/step)
    lat_index = int((LAT+45)/step)
    map_index = lon_index+lat_index*72+ascend_flag*1296
    return map_index

'''Give the index for bkg map for an array'''
def bkgmap_time(lon_arr,lat_arr,ascend_flag,bkgmap_flag):
    fnum = np.size(lon_arr)
    for ii in xrange(0,fnum):
        bkgmap_flag[ii] = bkgmap_index(lon_arr[ii],lat_arr[ii],ascend_flag[ii],5.)

'''Give the data points for specific 5x5 size'''
def bkgmap_num_cal(bkgmap_flag):
    dval = bkgmap_flag[1:(len(bkgmap_flag))] - bkgmap_flag[0:(len(bkgmap_flag)-1)]
    fnum = np.size((np.where(np.abs(dval) > 0))) + 1
    #print("5x5 interval number:",fnum)
    return fnum

'''Give the start and stop indices for specific 5x5 size'''
def bkgmap_interval_index(bkgmap_flag,bkgmap_start_index,bkgmap_stop_index):
    fnum = np.size(bkgmap_flag)
    inum = np.size(bkgmap_start_index)
    dval = bkgmap_flag[1:fnum] - bkgmap_flag[0:(fnum-1)]
    nonzero_in = np.where(np.abs(dval) > 0)
    nonzero_num = np.size(nonzero_in)
    #print( 'Non zeros index',nonzero_in[0:(fnum)] )
    bkgmap_start_index[1:(inum)] = nonzero_in + np.ones(nonzero_num,dtype=np.int)
    bkgmap_start_index[0] = 0
    #print( np.size(bkgmap_stop_index[0:(inum-2)]),np.size(nonzero_in) )
    bkgmap_stop_index[0:(inum-1)] =nonzero_in + np.zeros(nonzero_num,dtype=np.int)
    bkgmap_stop_index[inum-1] = fnum-1

def write_bkgspec(fname,channel,counts,spec_err,expo,hdr_ext):
    spec_qua= np.zeros(256)
    spec_grp= np.ones(256)
    spec_col1 = pf.Column(name='CHANNEL', format='J', array=channel)
    #spec_col2 = pf.Column(name='COUNTS', format='D', array=counts)
    spec_col2 = pf.Column(name='RATE', format='D',unit='counts/s', array=counts/expo)
    spec_col3 = pf.Column(name='STAT_ERR', format='D', array=spec_err/expo)
    spec_col4 = pf.Column(name='QUALITY', format='I', array=spec_qua)
    spec_col5 = pf.Column(name='GROUPING', format='I', array=spec_grp)
    cols = pf.ColDefs([spec_col1, spec_col2, spec_col3, spec_col4, spec_col5])
    hdr = pf.Header()
    hdr['EXTNAME']  = "SPECTRUM"
    hdr['PHAVERSN'] = "1992a"
    hdr['HDUCLASS'] = "OGIP"
    hdr['HDUCLAS1'] = "SPECTRUM"
    hdr['HDUCLAS2'] = "TOTAL"
    hdr['HDUCLAS3'] = "RATE"
    hdr['HDUCLAS4'] = "TYPE:I"
    hdr['HDUVERS1'] = "1.2.1"
    hdr['CHANTYPE'] = "PI"
    hdr['DETCHANS'] = hdr_ext['DETCHANS']
    hdr['TELESCOP'] = 'HXMT'
    hdr['INSTRUME'] = 'HE'

    hdr["OBS_MODE"] = hdr_ext['OBS_MODE']
    hdr["DATE-OBS"] = hdr_ext['DATE-OBS']
    hdr["DATE-END"] = hdr_ext['DATE-END']
    hdr["OBJECT"]   = hdr_ext['OBJECT']
    hdr["TSTART"]   = hdr_ext['TSTART']
    hdr["TSTOP"]    = hdr_ext['TSTOP']
    hdr["RA_OBJ"]   = hdr_ext['RA_OBJ']
    hdr["DEC_OBJ"]  = hdr_ext['DEC_OBJ']

    hdr['CORRFILE'] = 'None'
    hdr['CORRSCAL'] = 1.0

    hdr['BACKFILE'] = 'NONE'
    hdr['BACKSCAL'] = 1.0

    hdr['RESPFILE'] = 'NONE'
    hdr['ANCRFILE'] = 'NONE'
    hdr['FILTER']   = 'NONE'

    hdr['AREASCAL'] = 1.0
    hdr['EXPOSURE'] = expo
    #hdr['EXPOSURE'] = 1.0
    hdr['LIVETIME'] = expo
    hdr['DEADC']    = 1.0
    hdr['STATERR']  = True
    hdr['SYSERR']   = False
    hdr['POISSERR'] = False
    hdr['GROUPING'] = 0
    hdr['QUALITY']  = 0
    hdr['MJDREFI']  = hdr_ext['MJDREFI']
    hdr['MJDREFF']  = hdr_ext['MJDREFF']
    history_num = np.size(hdr_ext['HISTORY'])
    for ii in range(0,history_num):
        hdr['HISTORY']  = hdr_ext['HISTORY'][ii]
    hdu = pf.BinTableHDU.from_columns(cols,header=hdr)
    #hdu.writeto(fname)
    #hdu.writeto(fname,overwrite=True)
    try:
        hdu.writeto(fname,overwrite=True)
    except overwriteError:
        hdu.writeto(fname)

def write_lcurve(fname,time,counts,error,hdr_ext):
    lc_frac= np.ones(np.size(time))
    lc_col1 = pf.Column(name='Time', format='D', unit='s',array=time)
    lc_col2 = pf.Column(name='Rate', format='D', unit='counts/s', array=counts)
    lc_col3 = pf.Column(name='Error', format='D', array=error)
    lc_col4 = pf.Column(name='FRACEXP', format='D', array=lc_frac)
    cols = pf.ColDefs([lc_col1, lc_col2, lc_col3, lc_col4])
    hdr = pf.Header()
    hdr['EXTNAME']  = "RATE"
    hdr['PHAVERSN'] = "1992a"
    hdr['HDUCLASS'] = "OGIP"
    hdr['HDUCLAS1'] = "LIGHTCURVE"
    hdr['HDUCLAS2'] = "ALL"
    hdr['HDUCLAS3'] = "COUNT"
    hdr['HDUCLAS4'] = "TYPE:I"
    hdr['HDUVERS1'] = "1.1.0"
    hdr['CHANTYPE'] = "PI"
    hdr['TELESCOP'] = hdr_ext['TELESCOP']
    hdr['INSTRUME'] = hdr_ext['INSTRUME']
    hdr['TIMEUNIT'] = hdr_ext['TIMEUNIT']
    hdr['MJDREFI']  = hdr_ext['MJDREFI']
    hdr['MJDREFF']  = hdr_ext['MJDREFF']
    hdr["OBS_MODE"] = hdr_ext['OBS_MODE']
    hdr["DATE-OBS"] = hdr_ext['DATE-OBS']
    hdr["DATE-END"] = hdr_ext['DATE-END']
    hdr["OBJECT"]   = hdr_ext['OBJECT']
    hdr["TSTART"]   = hdr_ext['TSTART']
    hdr["TSTOP"]    = hdr_ext['TSTOP']
    hdr["RA_OBJ"]   = hdr_ext['RA_OBJ']
    hdr["DEC_OBJ"]  = hdr_ext['DEC_OBJ']
    hdr['TIMEZERO']  = hdr_ext['TIMEZERO']
    hdr['TIMEDEL']  = hdr_ext['TIMEDEL']
    history_num = np.size(hdr_ext['HISTORY'])
    for ii in range(0,history_num):
        hdr['HISTORY']  = hdr_ext['HISTORY'][ii]
    hdu = pf.BinTableHDU.from_columns(cols,header=hdr)
    #hdu.writeto(fname)
    #hdu.writeto(fname,overwrite=True)
    try:
        hdu.writeto(fname,overwrite=True)
    except overwriteError:
        hdu.writeto(fname)
# Gauss filter for pulsar profile

def GaussFilter(x,y,a):
    data_num = np.size(x)
    p        = np.ones(np.size(y))
    W        = np.zeros(np.size(y))
    yfilter  = np.zeros(np.size(y))
    for ii in range(data_num):
        W[ii]=np.sum(p*np.exp(-(x-x[ii])**2/(2*a**2)))

    for ii in range(data_num):
        yfilter[ii]=np.sum(p*y*np.exp(-(x-x[ii])**2/(2*a**2)))/W[ii]
    return yfilter

def Read_GTI(gtifile):
    '''Read gti extesion'''
    try:
        hdulist = pf.open(gtifile)
        tb = hdulist[1].data
        START = tb.field(0)
        STOP = tb.field(1)
        hdulist.close()
        print("GTI START=",START)
        print("GTI STOP=",STOP)
        return START,STOP
    except IOError:
        print(gtifile," is empty!")
        sys.exit()


def Read_deadtime(dtname,START,STOP,bkgmap_start_time,bkgmap_stop_time):
    bkgmap_num = np.size(bkgmap_start_time)
    '''Read Dead-Time file'''
    try:
        dtimelist = pf.open(dtname)
        dtime_tab = dtimelist[1].data
        dtime_time = dtime_tab.field(0)
        dtime_cyc = dtime_tab.field(1)
        dtime_num = np.size(dtime_time)
        dtime_arr = np.zeros((dtime_num,18))
        dt_cyc = np.zeros(dtime_num)
        for ii in xrange(0,18):
            if ii<=5:
                dt_cyc = dtime_cyc[0:(dtime_num),0]
            if ii>=6 & ii <=11:
                dt_cyc = dtime_cyc[0:(dtime_num),1]
            if ii>=12:
                dt_cyc = dtime_cyc[0:(dtime_num),2]
            tmpdt = dtime_tab.field(ii+2)/dt_cyc
            dtime_arr[0:(dtime_num),ii] = tmpdt
        dtimelist.close()
    except IOError:
        print(dtname," is empty!")
        sys.exit()

    #dtf = 1-np.mean(dtime_arr[0:(dtime_num),16])

    dtime_num = np.size(dtime_time)
    tflag = np.zeros(dtime_num,dtype='int')
    tindex = np.zeros(dtime_num,dtype='int')
    time_gtiflag(dtime_time,START,STOP,tflag)
    dtime_num = flag_selection(tflag,tindex)
    dtime_time=dtime_time[tindex[0:(dtime_num)]]
    dtime_cyc=dtime_cyc[tindex[0:(dtime_num)],0:3]
    dtime_arr=dtime_arr[tindex[0:(dtime_num)],0:18]

    bkgmap_dtc = np.zeros((18,bkgmap_num))
    bkgmap_dtt = np.zeros((18,bkgmap_num))

    #print(np.size(dtime_time))
    #print(dtime_num)

    '''Cal the dead time correction for every time interval'''
    for ii in xrange(0,18):
        for jj in xrange(0,bkgmap_num):
            t0 = bkgmap_start_time[jj]
            t1 = bkgmap_stop_time[jj]
            tflag2 = np.zeros(dtime_num,dtype='int')
            tindex2 = np.zeros(dtime_num,dtype='int')
            #print(np.size(dtime_time),dtime_num,dtime_time[0]-t0,dtime_time[dtime_num-1]-t1)
            time_gtiflag(dtime_time,t0,t1,tflag2)
            tmpdt_num = flag_selection(tflag2,tindex2)
            tmpdeadtime = dtime_arr[tindex2[0:tmpdt_num],ii]
            bkgmap_dtc[ii,jj] = np.sum(tmpdeadtime)/np.size(tmpdeadtime)
            bkgmap_dtt[ii,jj] = np.sum(tmpdeadtime)
            if np.size(tmpdeadtime)==0:
                bkgmap_dtc[ii,jj] = 0

    return bkgmap_dtc,bkgmap_dtt



'''Read Decay component file'''
'''And calculate the correction for decay component'''

def Read_decay(decname,tt_obs0):

    #decname = REFPATH + 'HE_decay.fits'

    try:
        dclist = pf.open(decname)
        dc_tab = dclist[1].data
        dc_index = dc_tab.field(0)
        dc_time = dc_tab.field(1)
        dc_spec = dc_tab.field(2)
        dclist.close()
    except IOError:
        print(decname," is empty!")
        sys.exit()

    dc_num = np.size(dc_index)
    
    dc_corr = np.zeros(4608);
    for ii in xrange(0,4608):
        tmpdc_spec = dc_spec[0:dc_num,ii]
        dc_corr[ii] = np.interp(tt_obs0,dc_time,tmpdc_spec)
    return dc_corr

'''Read DETID==17 event data (Blind detecter)'''

def Read_blind_detector(evtfilename):

    try:
        evt_list  = pf.open(evtfilename)
        evt_tab   = evt_list[1].data
        evt_time  = evt_tab.field(0)
        evt_detid = evt_tab.field(1)
        evt_cha   = evt_tab.field(2)
        evt_type  = evt_tab.field(5)
        evt_hdr   = evt_list[1].header
        evt_list.close()
    except IOError:
        print(evtfilename," is empty!")
        sys.exit()

    phy_index = np.where((evt_type == 0))
    evt_time_all  = evt_time[phy_index]
    evt_cha_all   = evt_cha[phy_index]
    evt_detid_all = evt_detid[phy_index]
    detid17_index = np.where((evt_detid == 16) & (evt_type == 0))
    evt_time  = evt_time[detid17_index]
    evt_cha   = evt_cha[detid17_index]

    return evt_time,evt_cha,evt_time_all,evt_cha_all,evt_detid_all

'''Read coe correct channel range'''
def Read_mchran(cha_data):

    cha_data_st = np.zeros((18,6))
    cha_data_sp = np.zeros((18,6))

    for ii in xrange(0,18):
        for jj in xrange(0,6):
            if (jj==0):
                cha_data_st[ii,jj] = int(cha_data[ii,jj])
            if (jj>0):
                cha_data_st[ii,jj] = int(cha_data[ii,jj])+1

    for ii in xrange(0,18):
        for jj in xrange(1,7):
            cha_data_sp[ii,jj-1] = int(cha_data[ii,jj])
    return cha_data_st,cha_data_sp

'''Obtain the blind spectra for each 5x5 degrees'''


def Calculate_blind_spectra(evt_time,evt_cha,bkgmap_dtc,bkgmap_start_time,bkgmap_stop_time):
    bkgmap_num = np.size(bkgmap_start_time)
    cha_ran = np.linspace(0,256,257)
    channel = np.linspace(0,255,256)
    #bld_spec_all_expo=0
    #bld_spec_all_deadt=np.sum(bkgmap_dtt[16,0:bkgmap_num])
    #bld_spec_all=np.zeros(256)
    bld_spec_arr=np.zeros((bkgmap_num,259),dtype='float')
    bld_dc_arr=np.zeros((bkgmap_num,3),dtype='float')

    for jj in xrange(0,bkgmap_num):
        t0 = bkgmap_start_time[jj]
        t1 = bkgmap_stop_time[jj]
        ttindex = np.where((evt_time >= t0) & (evt_time < t1+1))
        tmpcha = evt_cha[ttindex]
        #print( 'Exposure: ', t1-t0+1, ' Photon number: ',np.size(tmpcha))
        tmpspec,bins = np.histogram(tmpcha,bins=cha_ran,range=[0,255])
        #bld_spec_all_expo=bld_spec_all_expo+t1-t0+1
        #bld_spec_all=bld_spec_all+tmpspec
        bld_spec_arr[jj,0] = 16
        bld_spec_arr[jj,1] = (t0+t1)/2.0
        bld_spec_arr[jj,2] = t1-t0+1
        tmpdc = bkgmap_dtc[16,jj]
        bld_spec_arr[jj,3:259] = tmpspec/(1-tmpdc)
        '''
        if jj<54:
            continue
        plt.figure()
        plt.plot(channel,tmpspec)
        plt.show()
        print(tmpdc)
        print(tmpspec)
        print(tmpspec/(1-tmpdc))
        '''
    return bld_spec_arr

'''
Cal the spectrum for all detectors for correct the background lightcurve
'''

def Calculate_spectra_alldec(evt_time_all,evt_cha_all,bkgmap_dtc,bkgmap_start_time,bkgmap_stop_time):
    cha_ran = np.linspace(0,256,257)
    channel = np.linspace(0,255,256)

    bkgmap_num = np.size(bkgmap_start_time)
    spec_arr_all = np.zeros((18,256),dtype='float')
    dtcarr_all = np.zeros((18),dtype='float')

    for detider in xrange(0,18):
        for jj in xrange(0,bkgmap_num):
            t0 = bkgmap_start_time[jj]
            t1 = bkgmap_stop_time[jj]
            ttindex = np.where((evt_time_all >= t0) & (evt_time_all < t1+1) & (evt_detid_all==detider))
            tmpcha = evt_cha_all[ttindex]
            #print( 'Exposure: ', t1-t0+1, ' Photon number: ',np.size(tmpcha))
            tmpspec,bins = np.histogram(tmpcha,bins=cha_ran,range=[0,255])
            tmpdc = bkgmap_dtc[detider,jj]
            spec_arr_all[detider,0:256] = spec_arr_all[detider,0:256] + tmpspec
            dtcarr_all[detider]   = dtcarr_all[detider] + tmpdc

    return spec_arr_all,dtcarr_all

'''Read background model'''
def Read_background_model(srcmapname,srcerrmapname):

    try:
        srclist = pf.open(srcmapname)
        srcmap_tab = srclist[1].data
        srcmap_IN  = srcmap_tab.field(0)
        srcmap_LON = srcmap_tab.field(1)
        srcmap_LAT = srcmap_tab.field(2)
        srcmap_BKG = srcmap_tab.field(3)
        srclist.close()
        #print(srcmap_IN)
        #print(srcmap_LON)
        #print(srcmap_LAT)
    except IOError:
        print(srcmapname," is missing.")
        sys.exit()

    '''Read errors background model'''
    #srcmapname=REFPATH+'HE_bkgmaperr.fits'
    srcmapname=REFPATH+'HE_bkgmaperr_v20190711.fits'
    try:
        srclist2 = pf.open(srcmapname)
        srcmer_tab = srclist2[1].data
        srcmer_IN  = srcmer_tab.field(0)
        srcmer_LON = srcmer_tab.field(1)
        srcmer_LAT = srcmer_tab.field(2)
        srcmer_BKG = srcmer_tab.field(3)
        srclist2.close()
    except IOError:
        print(srcmapname," is missing.")
        sys.exit()
    return srcmap_BKG,srcmer_BKG

'''Read EHK file and position devision by 5 x 5'''
def Read_EHKfile(ehkname,START,STOP):
    try:
        ehklist = pf.open(ehkname)
        ehk_tab = ehklist[1].data
        ehk_time= ehk_tab.field(0)
        ehk_LON = ehk_tab.field(8)
        ehk_LAT = ehk_tab.field(9)
        ehklist.close()
    except IOError:
        print(ehkname," is empty!")
        sys.exit()
    ehk_num = np.size(ehk_time)

    #Ascend: 0
    #Descend: 1
    ehk_diff = ehk_LAT[1:(ehk_num)] - ehk_LAT[0:(ehk_num-1)] 
    ehk_ascend_flag = np.zeros(ehk_num);
    ehk_ascend_flag[np.where(ehk_diff >= 0)] = 0
    ehk_ascend_flag[np.where(ehk_diff < 0)]  = 1
    ehk_ascend_flag[ehk_num-1]  = ehk_ascend_flag[ehk_num-2]


    tflag = np.zeros(ehk_num,dtype='int')
    tindex = np.zeros(ehk_num,dtype='int')
    time_gtiflag(ehk_time,START,STOP,tflag)
    ehk_num = flag_selection(tflag,tindex)
    ehk_time = ehk_time[tindex[0:(ehk_num)]]
    ehk_LON = ehk_LON[tindex[0:(ehk_num)]]
    ehk_LAT = ehk_LAT[tindex[0:(ehk_num)]]
    ehk_ascend_flag= ehk_ascend_flag[tindex[0:(ehk_num)]]

    bkgmap_flag = np.zeros(ehk_num,dtype='int')
    bkgmap_time(ehk_LON,ehk_LAT,ehk_ascend_flag,bkgmap_flag)
    bkgmap_flag = np.array(bkgmap_flag)
    bkgmap_num = bkgmap_num_cal(bkgmap_flag)

    bkgmap_start_index = np.zeros(bkgmap_num,dtype='int')
    bkgmap_stop_index = np.zeros(bkgmap_num,dtype='int')
    bkgmap_interval_index(bkgmap_flag,bkgmap_start_index,bkgmap_stop_index)
    bkgmap_arr_num  = bkgmap_stop_index - bkgmap_start_index + 1
    bkgmap_arr_expo = bkgmap_stop_index - bkgmap_start_index + 1

    bkgmap_start_time = ehk_time[bkgmap_start_index]
    bkgmap_stop_time  = ehk_time[bkgmap_stop_index]
    bkgmap_flag_uniq  = bkgmap_flag[bkgmap_start_index]
    return bkgmap_start_time,bkgmap_stop_time,bkgmap_flag_uniq,bkgmap_num,bkgmap_arr_expo


'''Read EHK file and position devision by 5 x 5'''
def Read_EHKfile_lc(ehkname,START,STOP):
    try:
        ehklist = pf.open(ehkname)
        ehk_tab = ehklist[1].data
        ehk_time= ehk_tab.field(0)
        ehk_LON = ehk_tab.field(8)
        ehk_LAT = ehk_tab.field(9)
        ehklist.close()
    except IOError:
        print(ehkname," is empty!")
        sys.exit()
    ehk_num = np.size(ehk_time)

    #Ascend: 0
    #Descend: 1
    ehk_diff = ehk_LAT[1:(ehk_num)] - ehk_LAT[0:(ehk_num-1)] 
    ehk_ascend_flag = np.zeros(ehk_num);
    ehk_ascend_flag[np.where(ehk_diff >= 0)] = 0
    ehk_ascend_flag[np.where(ehk_diff < 0)]  = 1
    ehk_ascend_flag[ehk_num-1]  = ehk_ascend_flag[ehk_num-2]

    tflag = np.zeros(ehk_num,dtype='int')
    tindex = np.zeros(ehk_num,dtype='int')
    time_gtiflag(ehk_time,START,STOP,tflag)
    ehk_num = flag_selection(tflag,tindex)
    ehk_time = ehk_time[tindex[0:(ehk_num)]]
    ehk_LON = ehk_LON[tindex[0:(ehk_num)]]
    ehk_LAT = ehk_LAT[tindex[0:(ehk_num)]]
    ehk_ascend_flag= ehk_ascend_flag[tindex[0:(ehk_num)]]

    delta_T = 16
    Tmin = np.min(ehk_time)
    Tmax = np.max(ehk_time)
    Tnum = int((Tmax-Tmin)/delta_T)
    T_arr_low = Tmin + np.arange(0,Tnum,1,dtype=np.float64)*delta_T
    T_arr_up  = Tmin + (np.arange(0,Tnum,1,dtype=np.float64)+1)*delta_T
    T_arr = (T_arr_low+T_arr_up)/2.0
    ehk_time2 = np.zeros(Tnum,dtype=float)
    ehk_LON2 = np.zeros(Tnum,dtype=float)
    ehk_LAT2 = np.zeros(Tnum,dtype=float)
    ehk_ascend_flag2=np.zeros(Tnum,dtype=float)
    expo_flag2=np.zeros(Tnum,dtype=int)
    bkgmap_arr_expo2=np.zeros(Tnum,dtype=int)
    for ii in range(0,Tnum):
        tmpt0 = T_arr_low[ii]
        tmpt1 = T_arr_up[ii]
        tmpin = np.where((ehk_time>=tmpt0)&(ehk_time<tmpt1))
        if (np.size(tmpin)==0):
            expo_flag2[ii] = 1
            continue
        ehk_time2[ii] = (tmpt0+tmpt1)/2.0
        ehk_LON2[ii] = np.mean(ehk_LON[tmpin])
        ehk_LAT2[ii] = np.mean(ehk_LAT[tmpin])
        bkgmap_arr_expo2[ii] = np.size(tmpin)
        if (bkgmap_arr_expo2[ii]<delta_T/2.):
            expo_flag2[ii] = 1

    tmpindex = np.where(expo_flag2 == 0)
    ehk_time2 = ehk_time2[tmpindex]
    ehk_LON2 = ehk_LON2[tmpindex]
    ehk_LAT2 = ehk_LAT2[tmpindex]
    bkgmap_arr_expo2=bkgmap_arr_expo2[tmpindex]
    Tnum = np.size(ehk_time2)

    ehk_diff2 = ehk_LAT2[1:(Tnum)] - ehk_LAT2[0:(Tnum-1)] 
    ehk_ascend_flag2[np.where(ehk_diff2 >= 0)] = 0
    ehk_ascend_flag2[np.where(ehk_diff2 < 0)]  = 1
    ehk_ascend_flag2[Tnum-1]  = ehk_ascend_flag2[Tnum-2]

    bkgmap_flag = np.zeros(Tnum,dtype='int')
    bkgmap_time(ehk_LON2,ehk_LAT2,ehk_ascend_flag2,bkgmap_flag)
    bkgmap_flag = np.array(bkgmap_flag)

    bkgmap_start_time = ehk_time2-delta_T/2.0
    bkgmap_stop_time  = ehk_time2+delta_T/2.0
    bkgmap_flag_uniq  = bkgmap_flag
    bkgmap_num = np.size(bkgmap_start_time)
    bkgmap_arr_expo = bkgmap_arr_expo2
    print(bkgmap_start_time-bkgmap_start_time[0])
    print(bkgmap_arr_expo2)
    return bkgmap_start_time,bkgmap_stop_time,bkgmap_flag_uniq,bkgmap_num,bkgmap_arr_expo



'''Cal the spectra of all detectors from map'''

def Calculate_background(bkgmap_flag_uniq,bkgmap_start_time,bkgmap_stop_time,bkgmap_arr_expo,srcmap_BKG,srcmer_BKG,dc_corr,bld_spec_arr):

    bkgspec_all_bkgmap = np.zeros((18*bkgmap_num,259),dtype='float')
    bkgspec_err_bkgmap = np.zeros((18*bkgmap_num,259),dtype='float')

    for detid in xrange(0,18):
        tin1 = detid*256
        tin2 = (detid+1)*256
        for ii in xrange(0,bkgmap_num):
            tmpindex = bkgmap_flag_uniq[ii]
            tmpspec  = srcmap_BKG[tmpindex,0:4608] - dc_corr
            bkgspec_all_bkgmap[detid*bkgmap_num+ii,0] = detid
            bkgspec_all_bkgmap[detid*bkgmap_num+ii,1] = (bkgmap_start_time[ii]+bkgmap_start_time[ii])/2.
            bkgspec_all_bkgmap[detid*bkgmap_num+ii,2] = bkgmap_arr_expo[ii]
            bkgspec_all_bkgmap[detid*bkgmap_num+ii,3:259] = tmpspec[tin1:tin2]*bkgmap_arr_expo[ii]
            ttspec = tmpspec[tin1:tin2]
            inin = np.where(ttspec > 950)
            ttspec[inin] = ttspec[inin] - 999
            bkgspec_all_bkgmap[detid*bkgmap_num+ii,3:259] = ttspec*bkgmap_arr_expo[ii]
            tmpspecerr  = srcmer_BKG[tmpindex,0:4608]
            bkgspec_err_bkgmap[detid*bkgmap_num+ii,0] = detid
            bkgspec_err_bkgmap[detid*bkgmap_num+ii,1] = (bkgmap_start_time[ii]+bkgmap_start_time[ii])/2.
            bkgspec_err_bkgmap[detid*bkgmap_num+ii,2] = bkgmap_arr_expo[ii]
            bkgspec_err_bkgmap[detid*bkgmap_num+ii,3:259] = tmpspecerr[tin1:tin2]*bkgmap_arr_expo[ii]
            #print("ii===",ii,tmpindex,"tmpspec==",np.sum(tmpspec[0:256]),np.sum(dc_corr))
            #print("ii===",ii,tmpindex,"tmpspec==",np.sum(tmpspec[0:256]),np.sum(dc_corr))

    '''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
    '''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
    '''+++++++++++++++++++++++++++++++++++++++Calculate spectrum fro BLD'''
    '''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
    '''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
    #spec_bldmod=np.zeros(256)
    #spec_bldmod_expo=0
    spec_ch = np.linspace(0,255,256)
    tmpexpo_arr = bkgspec_all_bkgmap[16*bkgmap_num:((16+1)*bkgmap_num),2]
    tmpspec_arr = bkgspec_all_bkgmap[16*bkgmap_num:((16+1)*bkgmap_num),3:259]
    #spec_bldmod= np.sum(tmpspec_arr[valid_flag_in],axis=0)
    #spec_bldmod_expo=np.sum(tmpexpo_arr[valid_flag_in])
    spec_bldmod_arr= bkgspec_all_bkgmap[16*bkgmap_num:((16+1)*bkgmap_num),0:259]

    #print(bld_spec_all_expo,spec_bldmod_expo)
    #print(np.sum(bld_spec_all),np.sum(spec_bldmod))

    #plt.figure()
    #plt.plot(channel,bld_spec_all/(1-bld_spec_all_deadc))
    #plt.plot(channel,bld_spec_all)
    #plt.plot(channel,rr)
    #plt.show()
    #plt.figure()
    #plt.plot(channel,spec_bldmod-(bld_spec_all/(1-bld_spec_all_deadc)))
    #plt.show()

    bld_corarr = np.zeros(6)

    bld_corarr = np.zeros((bkgmap_num,259),dtype='float')
    bld_corerr = np.zeros((bkgmap_num,259),dtype='float')

    bkgspec_all_bkgmap_cor = np.zeros((18*bkgmap_num,259),dtype='float')
    bkgspec_err_bkgmap_cor = np.zeros((18*bkgmap_num,259),dtype='float')
    bkgspec_bkgmap_cor_flag = np.zeros(bkgmap_num,dtype='float') + 1

    tmpr = np.zeros(256,dtype='float')
    for ii in xrange(0,bkgmap_num):
        bld_corarr[ii,0] = bld_spec_arr[ii,0]
        bld_corarr[ii,1] = bld_spec_arr[ii,1]
        bld_corarr[ii,2] = bld_spec_arr[ii,2]
        tmpr = bld_spec_arr[ii,3:259] / spec_bldmod_arr[ii,3:259]
        #print(np.sum(tmpr))
        bld_corarr[ii,3:259] = tmpr
        bld_corerr[ii,0] = bld_spec_arr[ii,0]
        bld_corerr[ii,1] = bld_spec_arr[ii,1]
        bld_corerr[ii,2] = bld_spec_arr[ii,2]
        bld_corerr[ii,3:259]=np.sqrt(bld_spec_arr[ii,3:259]) / spec_bldmod_arr[ii,3:259]
        if ((np.sum(tmpr)<=0) | (bld_spec_arr[ii,2]<=1.)):
            bkgspec_bkgmap_cor_flag[ii] = 0
        '''plt.figure()
        plt.plot(bld_spec_arr[ii,3:259])
        plt.plot(spec_bldmod_arr[ii,3:259])
        plt.plot(tmpr)
        plt.show()
        print(bld_spec_arr[ii,3:259])'''

    for detid in xrange(0,18):
        tin1 = detid*256
        tin2 = (detid+1)*256
        sindex = np.array(chrelation_bld[0:256,detid],dtype='int')
        for ii in xrange(0,bkgmap_num):
            tmpindex = detid*bkgmap_num+ii

            tmpr = bld_corarr[ii,3:259]
            tmpcorrerr  = bld_corerr[ii,3:259]

            tmpspec  = bkgspec_all_bkgmap[tmpindex,3:259]
            tmpspecerr  = bkgspec_err_bkgmap[tmpindex,3:259]

            tmpspec_cor = tmpspec*tmpr[sindex]
            tmpspecerr_cor = np.sqrt( (tmpr[sindex]*tmpspecerr)**2 + (tmpcorrerr[sindex]*tmpspec)**2)
            #plt.figure()
            #plt.plot(tmpspec)
            #plt.plot(tmpspec_cor)
            #plt.show()
            bkgspec_all_bkgmap_cor[tmpindex,0] = bkgspec_all_bkgmap[tmpindex,0]
            bkgspec_all_bkgmap_cor[tmpindex,1] = bkgspec_all_bkgmap[tmpindex,1]
            bkgspec_all_bkgmap_cor[tmpindex,2] = bkgspec_all_bkgmap[tmpindex,2]
            bkgspec_all_bkgmap_cor[tmpindex,3:259] = tmpspec_cor
            bkgspec_all_bkgmap_cor[tmpindex,12:24] = 0
            #print(np.sum(tmpr))
            #if (ii>= 0) & (detid==0):
                #if (np.sum(tmpspec_cor)==0):
                #    bkgspec_bkgmap_cor_flag[tmpindex] = 1
                #print('det ',detid,' seq ', ii, ' expo=',bkgspec_all_bkgmap[tmpindex,2], ' cnt= ',np.sum(tmpspec),np.sum(tmpspec_cor))
                #plt.figure()
                #plt.plot(tmpspec_cor)
                #plt.plot(tmpspec)
                #plt.plot(tmpr[sindex])
                #plt.show()
            bkgspec_err_bkgmap_cor[tmpindex,0] = bkgspec_err_bkgmap[tmpindex,0]
            bkgspec_err_bkgmap_cor[tmpindex,1] = bkgspec_err_bkgmap[tmpindex,1]
            bkgspec_err_bkgmap_cor[tmpindex,2] = bkgspec_err_bkgmap[tmpindex,2]
            bkgspec_err_bkgmap_cor[tmpindex,3:259] = tmpspecerr_cor
            bkgspec_err_bkgmap_cor[tmpindex,12:24] = 1
            #print("ii===",ii,tmpindex,"tmpspec==",np.sum(tmpspec[0:256]),np.sum(dc_corr))
            #print("ii===",ii,tmpindex,"tmpspec==",np.sum(tmpspec[0:256]),np.sum(dc_corr))
    return bkgspec_all_bkgmap_cor,bkgspec_err_bkgmap_cor,bkgspec_bkgmap_cor_flag


'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++  Calculate spectrum    +++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''

if sp_lc_select == 'spec':
    '''Read GTI file'''
    [START,STOP] = Read_GTI(gtifile)

    '''Read EHK file and position devision by 5 x 5'''
    bkgmap_start_time,bkgmap_stop_time,bkgmap_flag_uniq,bkgmap_num,bkgmap_arr_expo=Read_EHKfile(ehkname,START,STOP)

    '''Read dead time'''
    bkgmap_dtc,bkgmap_dtt = Read_deadtime(dtname,START,STOP,bkgmap_start_time,bkgmap_stop_time)

    '''Calculate the decay component'''
    decname = REFPATH + 'HE_decay_v20190711.fits'
    tt_obs0 = np.mean(bkgmap_start_time)
    dc_corr = Read_decay(decname,tt_obs0)

    '''Read blind detector events'''
    evt_time,evt_cha,evt_time_all,evt_cha_all,evt_detid_all = Read_blind_detector(evtfilename)
    cha_data = np.loadtxt(REFPATH+'mchran.txt')
    cha_data_st,cha_data_sp = Read_mchran(cha_data)

    '''Read channel relation between blind detector and the rest detectors'''
    chrelation_bld = np.loadtxt(REFPATH+'pha_channel_correspond.txt')

    '''Obtain the blind spectra for each 5x5 degrees'''
    bld_spec_arr = Calculate_blind_spectra(evt_time,evt_cha,bkgmap_dtc,bkgmap_start_time,bkgmap_stop_time)

    '''
    Cal the spectrum for all detectors for correct the background lightcurve
    '''

    spec_arr_all,dtcarr_all = Calculate_spectra_alldec(evt_time_all,evt_cha_all,bkgmap_dtc,bkgmap_start_time,bkgmap_stop_time)

    '''Read background model'''
    srcmapname=REFPATH+'HE_bkgmap_v20190711.fits'
    srcerrmapname=REFPATH+'HE_bkgmaperr_v20190711.fits'

    srcmap_BKG,srcmer_BKG=Read_background_model(srcmapname,srcerrmapname)

    '''Cal the spectra of blind detectors from map'''

    bkgspec_bld_bkgmap = np.zeros((bkgmap_num,256),dtype='float')
    valid_flag = np.zeros(bkgmap_num)+1
    for ii in xrange(0,bkgmap_num):
        tmpindex = bkgmap_flag_uniq[ii]
        tmpspec  = srcmap_BKG[tmpindex,0:4608]
        tin1 = 16*256
        tin2 = 17*256
        bkgspec_bld_bkgmap[ii,0:256] = tmpspec[tin1:tin2]
        if(np.mean(tmpspec[tin1:tin2])>900):
            valid_flag[ii] = 0


    #plt.figure()
    #plt.plot(dc_corr)
    #plt.show()
    '''Cal the spectra of all detectors from map'''

    bkgspec_all_bkgmap_cor,bkgspec_err_bkgmap_cor,bkgspec_bkgmap_cor_flag = Calculate_background(bkgmap_flag_uniq,bkgmap_start_time,bkgmap_stop_time,bkgmap_arr_expo,srcmap_BKG,srcmer_BKG,dc_corr,bld_spec_arr)

    ffgg = valid_flag * bkgspec_bkgmap_cor_flag
    valid_flag_in = np.where(ffgg == 1)
    #cor_flag_in = np.where(bkgspec_bkgmap_cor_flag == 1)

    '''Estimate background spectra'''
    spec_id   = np.zeros(18)-1
    src_name=[]
    #print(sl_name)
    sf = open(sl_name)
    cnter = 0
    for sline in sf:
        #print(sline)
        src_name.append(sline)
        tmppos = sline.find('\n')
        if (len(sline[0:tmppos]) == 0):
            print("Input file name error:")
            sys.exit()
        tspec_list    = pf.open(sline[0:tmppos])
        tspec_tab     = tspec_list[3].data
        spec_id[cnter]   = tspec_tab.field(0)
        tspec_list.close()
        cnter = cnter+1
    sf.close()
    if len(src_name) != 18 | cnter != 18:
        print("The input file number is not 18!")
        sys.exit()

    print(src_name)
    print("Calculate background spectra now.")
    spec_ch = np.linspace(0,255,256)
    specnum = np.size(src_name)
    for ii in xrange(0,specnum):
        iidd = int(spec_id[ii])
        tmpstr = src_name[ii]
        tmppos = tmpstr.find('\n')
        if (len(tmpstr[0:tmppos]) == 0):
            print("Input file name error:")
            sys.exit()
        spec_list    = pf.open(tmpstr[0:tmppos])
        spec_tab     = spec_list[1].data
        spec_channel = spec_tab.field(0)
        spec_counts  = spec_tab.field(1)
        spec_hdr     = spec_list[1].header
        spec_tab     = spec_list[3].data
        #spec_id      = spec_tab.field(0)
        spec_list.close()
        print("For detector:", iidd, " src file name: ", tmpstr[0:tmppos])
        tmpexpo_arr = bkgspec_all_bkgmap_cor[iidd*bkgmap_num:((iidd+1)*bkgmap_num),2]
        tmpspec_arr = bkgspec_all_bkgmap_cor[iidd*bkgmap_num:((iidd+1)*bkgmap_num),3:259]
        tmpspec_err = bkgspec_err_bkgmap_cor[iidd*bkgmap_num:((iidd+1)*bkgmap_num),3:259]
        tmpexpo = 0
        spec_cnt = np.zeros(256,dtype='float')
        spec_err = np.zeros(256,dtype='float')
        tmpexpo=np.sum(tmpexpo_arr[valid_flag_in])
        spec_cnt= np.sum(tmpspec_arr[valid_flag_in],axis=0)
        spec_err= np.sqrt(np.sum(tmpspec_err[valid_flag_in]*tmpspec_err[valid_flag_in],axis=0))
       
        tmpspec0 = spec_counts/spec_hdr['exposure']
        tmpspec1 = spec_cnt/tmpexpo
        #cnt0 = (np.sum(tmpspec0[220:256]) + np.sum(tmpspec0[0:10]))
        #cnt1 = (np.sum(tmpspec1[220:256]) + np.sum(tmpspec1[0:10]))
        ttspec0 = tmpspec0[220:256]
        ttspec1 = tmpspec1[220:256]
        cnt0 = (np.sum(tmpspec0[220:256]))
        cnt1 = (np.sum(tmpspec1[220:256]))
        if iidd==14:
            cnt0 = (np.sum(tmpspec0[220:253]))
            cnt1 = (np.sum(tmpspec1[220:253]))

        rr0  = cnt0/cnt1

        #rr = (spec_counts/spec_hdr['exposure'])/(spec_cnt/tmpexpo)
        #rr0=(np.sum(rr[220:256]) + np.sum(rr[0:10]))/46.0
        #print(rr0,cnt0,cnt1,tmpexpo)
        outname = outnam + '_' + str(iidd)+'.pha'
        write_bkgspec(outname,spec_ch,spec_cnt*rr0,spec_err*rr0,tmpexpo,spec_hdr)

        if(ii >= 20):
            print(rr0)
            plt.figure()
            plt.plot(spec_counts/spec_hdr['exposure'],'C1')
            plt.plot(spec_cnt/tmpexpo,'C2')
            plt.plot(spec_cnt/tmpexpo*rr0,'C3')
            #plt.plot(spec_cnt-spec_bldmod)
            plt.show()

        #plt.figure()
        #print("Expo=",tmpexpo)
        #plt.plot(rr)
        #plt.show()

'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++  Calculate light curve    +++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''


if sp_lc_select == 'lc':

    '''Read GTI file'''
    [START,STOP] = Read_GTI(gtifile)

    '''Read EHK file and position devision by 5 x 5'''
    bkgmap_start_time,bkgmap_stop_time,bkgmap_flag_uniq,bkgmap_num,bkgmap_arr_expo=Read_EHKfile_lc(ehkname,START,STOP)

    '''Read dead time'''
    bkgmap_dtc,bkgmap_dtt = Read_deadtime(dtname,START,STOP,bkgmap_start_time,bkgmap_stop_time)

    '''Calculate the decay component'''
    decname = REFPATH + 'HE_decay_v20190711.fits'
    tt_obs0 = np.mean(bkgmap_start_time)
    dc_corr = Read_decay(decname,tt_obs0)

    '''Read blind detector events'''
    evt_time,evt_cha,evt_time_all,evt_cha_all,evt_detid_all = Read_blind_detector(evtfilename)
    cha_data = np.loadtxt(REFPATH+'mchran.txt')
    cha_data_st,cha_data_sp = Read_mchran(cha_data)

    '''Read channel relation between blind detector and the rest detectors'''
    chrelation_bld = np.loadtxt(REFPATH+'pha_channel_correspond.txt')

    '''Obtain the blind spectra for each 5x5 degrees'''
    bld_spec_arr = Calculate_blind_spectra(evt_time,evt_cha,bkgmap_dtc,bkgmap_start_time,bkgmap_stop_time)

    '''
    Cal the spectrum for all detectors for correct the background lightcurve
    '''

    spec_arr_all,dtcarr_all = Calculate_spectra_alldec(evt_time_all,evt_cha_all,bkgmap_dtc,bkgmap_start_time,bkgmap_stop_time)

    '''Read background model'''
    srcmapname=REFPATH+'HE_bkgmap_v20190711.fits'
    srcerrmapname=REFPATH+'HE_bkgmaperr_v20190711.fits'

    srcmap_BKG,srcmer_BKG=Read_background_model(srcmapname,srcerrmapname)

    '''Cal the spectra of blind detectors from map'''

    bkgspec_bld_bkgmap = np.zeros((bkgmap_num,256),dtype='float')
    valid_flag = np.zeros(bkgmap_num)+1
    for ii in xrange(0,bkgmap_num):
        tmpindex = bkgmap_flag_uniq[ii]
        tmpspec  = srcmap_BKG[tmpindex,0:4608]
        tin1 = 16*256
        tin2 = 17*256
        bkgspec_bld_bkgmap[ii,0:256] = tmpspec[tin1:tin2]
        if(np.mean(tmpspec[tin1:tin2])>900):
            valid_flag[ii] = 0


    #plt.figure()
    #plt.plot(dc_corr)
    #plt.show()
    '''Cal the spectra of all detectors from map'''

    bkgspec_all_bkgmap_cor,bkgspec_err_bkgmap_cor,bkgspec_bkgmap_cor_flag = Calculate_background(bkgmap_flag_uniq,bkgmap_start_time,bkgmap_stop_time,bkgmap_arr_expo,srcmap_BKG,srcmer_BKG,dc_corr,bld_spec_arr)

    ffgg = valid_flag * bkgspec_bkgmap_cor_flag
    valid_flag_in = np.where(ffgg == 1)
    #cor_flag_in = np.where(bkgspec_bkgmap_cor_flag == 1)


    '''Obtain the correction factor first'''
    coe_coe = np.zeros(18)+1
    spec_ch = np.linspace(0,255,256)
    specnum = 18
    for ii in xrange(0,specnum):
        iidd = ii
        tmpexpo_arr = bkgspec_all_bkgmap_cor[iidd*bkgmap_num:((iidd+1)*bkgmap_num),2]
        tmpspec_arr = bkgspec_all_bkgmap_cor[iidd*bkgmap_num:((iidd+1)*bkgmap_num),3:259]
        tmpspec_err = bkgspec_err_bkgmap_cor[iidd*bkgmap_num:((iidd+1)*bkgmap_num),3:259]
        tmpexpo = 0
        spec_cnt = np.zeros(256,dtype='float')
        spec_err = np.zeros(256,dtype='float')
        tmpexpo=np.sum(tmpexpo_arr)
        spec_cnt= np.sum(tmpspec_arr,axis=0)
        spec_err= np.sqrt(np.sum(tmpspec_err*tmpspec_err,axis=0))
        
        tmpspec0 = spec_arr_all[iidd,0:256]/(tmpexpo-dtcarr_all[ii])
        tmpspec1 = spec_cnt/tmpexpo
        cnt0 = (np.sum(tmpspec0[220:256]))
        cnt1 = (np.sum(tmpspec1[220:256]))
        if iidd==14:
            cnt0 = (np.sum(tmpspec0[220:253]))
            cnt1 = (np.sum(tmpspec1[220:253]))

        rr0  = cnt0/cnt1
        #print(rr0,cnt0,cnt1,tmpexpo)
        if cnt0 >0:
            coe_coe[ii] = rr0
        if(ii >= 20):
            plt.figure()
            plt.plot(spec_counts/spec_hdr['exposure'],'C1')
            plt.plot(spec_cnt/tmpexpo,'C2')
            plt.plot(spec_cnt/tmpexpo*rr0,'C3')
            #plt.plot(spec_cnt-spec_bldmod)
            plt.show()
    '''Cal bakcground lightcurve!'''
    print("Calculate background lightcurve now.")
    src_name=[]
    print(sl_name)
    sf = open(sl_name)
    for sline in sf:
        #print( sline )
        src_name.append(sline)
    sf.close()
    #print(src_name)
    tmpstr = src_name[0]
    tmppos = tmpstr.find('\n')
    if (len(tmpstr[0:tmppos]) == 0):
        print("Input file name error:")
        sys.exit()
    #print("For HE src file name: ", tmpstr[0:tmppos])
    lc_list    = pf.open(tmpstr[0:tmppos])
    lc_tab     = lc_list[1].data
    lc_time    = lc_tab.field(0)
    lc_counts  = lc_tab.field(1)
    lc_hdr     = lc_list[1].header
    lc_list.close()
    lc_num = np.size(lc_time)
    lc_bkg     = np.zeros(lc_num)
    lc_bkg_map = np.zeros(bkgmap_num)
    lc_bkg_err = np.zeros(bkgmap_num)
    lc_tim_map=bkgspec_all_bkgmap_cor[0:bkgmap_num,1]
    #if(chmin <0):
    #    print("Illigal minmum channel, which will be set to 0")
    #    chmin=0
    #if(chmax > 255):
    #    print("Illigal maximum channel, which will be set to 255")
    #    chmmax=255
    #
    for idid in xrange(0,18):
        if (idid == 16):
            continue
        #print("idid== ",idid)
        for jj in xrange(0,bkgmap_num):
            if((valid_flag[jj]==0)|(bkgspec_bkgmap_cor_flag[jj] == 0)):
                lc_bkg_map[jj] = -1e6
                continue
            bkgindex0 = idid*bkgmap_num + jj
            tmpchmin = chmin[idid] + 3
            tmpchmax = chmax[idid] + 3 + 1
            #print(tmpchmin,tmpchmax)
            tmpexpo_arr = bkgspec_all_bkgmap_cor[bkgindex0,2]
            tmpspec_arr = bkgspec_all_bkgmap_cor[bkgindex0,tmpchmin:tmpchmax]
            tmpspec_err = bkgspec_err_bkgmap_cor[bkgindex0,tmpchmin:tmpchmax]
            spec_cnt= np.sum(tmpspec_arr) * coe_coe[idid]
            spec_err= np.sqrt(np.sum(tmpspec_err*tmpspec_err))/np.sum(tmpexpo_arr) * coe_coe[idid]
            #print(spec_cnt/np.sum(tmpexpo_arr),tmpexpo_arr)
            lc_bkg_map[jj] = lc_bkg_map[jj] + spec_cnt/np.sum(tmpexpo_arr)
            lc_bkg_err[jj] = np.sqrt(lc_bkg_err[jj]*lc_bkg_err[jj] + spec_err*spec_err)
            #plt.figure()
            #plt.plot(bkgspec_all_bkgmap_cor[bkgindex0,tmpchmin:256])
            #plt.show()
        #plt.figure()
        #plt.plot(lc_tim_map-lc_tim_map[0],lc_bkg_map)
        #plt.show()
    
    #fbkg=interpolate.interp1d(lc_tim_map,lc_bkg_map,kind='cubic')
    #lc_bkg=fbkg(lc_time)
    #print(lc_bkg_map)
    norm_in = np.where(lc_bkg_map >= -1e4)
    #print(lc_bkg_map[norm_in])
    tmptime_yf = lc_tim_map[norm_in]
    tmprate_yf = lc_bkg_map[norm_in]
    tmpyf_bkg0 = GaussFilter(tmptime_yf,tmprate_yf,100)
    tmpyf_res0 = tmprate_yf - tmpyf_bkg0
    tmpyf_res1 = GaussFilter(tmptime_yf,tmpyf_res0,200)
    tmprate_yf2 = tmpyf_bkg0 + tmpyf_res1
    lc_bkg = np.interp(lc_time,tmptime_yf,tmprate_yf2)
    lc_err = np.interp(lc_time,tmptime_yf,lc_bkg_err[norm_in])
    #lc_bkg = np.interp(lc_time,lc_tim_map[norm_in],lc_bkg_map[norm_in])
    #lc_err = np.interp(lc_time,lc_tim_map[norm_in],lc_bkg_err[norm_in])
    print(lc_bkg_map)
    print(lc_tim_map-lc_tim_map[0])
    print(lc_time-lc_tim_map[0])
    #plt.figure()
    #plt.plot(lc_tim_map[norm_in]-lc_tim_map[0],lc_bkg_map[norm_in])
    #plt.plot(lc_time-lc_tim_map[0],lc_bkg,'r-s')
    #plt.show()
    print("The total light curve: ")
    outname = outnam + '_all.lc'
    write_lcurve(outname,lc_time,lc_bkg,lc_err,lc_hdr)
    ''' For 0-17 '''

    for idid in xrange(0,18):
        if (idid == 16):
            continue
        #print("idid== ",idid)
        for jj in xrange(0,bkgmap_num):
            if((valid_flag[jj]==0)|(bkgspec_bkgmap_cor_flag[jj] == 0)):
                lc_bkg_map[jj] = -1e6
                continue
            bkgindex0 = idid*bkgmap_num + jj
            tmpchmin = chmin[idid] + 3
            tmpchmax = chmax[idid] + 3 + 1
            #print(tmpchmin,tmpchmax)
            tmpexpo_arr = bkgspec_all_bkgmap_cor[bkgindex0,2]
            tmpspec_arr = bkgspec_all_bkgmap_cor[bkgindex0,tmpchmin:tmpchmax]
            tmpspec_err = bkgspec_err_bkgmap_cor[bkgindex0,tmpchmin:tmpchmax]
            spec_cnt= np.sum(tmpspec_arr) * coe_coe[idid]
            spec_err= np.sqrt(np.sum(tmpspec_err*tmpspec_err))/np.sum(tmpexpo_arr) * coe_coe[idid]
            #print(spec_cnt/np.sum(tmpexpo_arr),tmpexpo_arr)
            lc_bkg_map[jj] = spec_cnt/np.sum(tmpexpo_arr)
            lc_bkg_err[jj] = spec_err*spec_err
            #plt.figure()
            #plt.plot(bkgspec_all_bkgmap_cor[bkgindex0,tmpchmin:256])
            #plt.show()
        #plt.figure()
        #plt.plot(lc_tim_map-lc_tim_map[0],lc_bkg_map)
        #plt.show()
    
        norm_in = np.where(lc_bkg_map >= -1e4)
        tmptime_yf = lc_tim_map[norm_in]
        tmprate_yf = lc_bkg_map[norm_in]
        tmpyf_bkg0 = GaussFilter(tmptime_yf,tmprate_yf,100)
        tmpyf_res0 = tmprate_yf - tmpyf_bkg0
        tmpyf_res1 = GaussFilter(tmptime_yf,tmpyf_res0,200)
        tmprate_yf2 = tmpyf_bkg0 + tmpyf_res1
        lc_bkg = np.interp(lc_time,tmptime_yf,tmprate_yf2)
        lc_err = np.interp(lc_time,tmptime_yf,lc_bkg_err[norm_in])
        #print(lc_bkg_map)
        #print(lc_tim_map-lc_tim_map[0])
        #print(lc_time-lc_tim_map[0])
        #plt.figure()
        #plt.plot(lc_tim_map[norm_in]-lc_tim_map[0],lc_bkg_map[norm_in])
        #plt.plot(lc_time-lc_tim_map[0],lc_bkg,'r-s')
        #plt.show()
        print("The total light curve: ")
        outname = outnam + '_ID' + str(idid) +'.lc'
        write_lcurve(outname,lc_time,lc_bkg,lc_err,lc_hdr)


'''
    for idid in xrange(0,18):
        if (idid == 16):
            continue
        print("For detector:", idid)
        for ii in xrange(0,lc_num):
            for jj in xrange(0,bkgmap_num):
                bkgindex0 = idid*bkgmap_num + jj
                tt0 = bkgmap_start_time[jj]
                tt1 = bkgmap_stop_time[jj] + 1
                tmpchmin = chmin + 3
                tmpchmax = chmax + 3
                if (lc_time[ii]>=tt0) & (lc_time[ii]< tt1):
                    tmpexpo_arr = bkgspec_all_bkgmap[bkgindex0,2]
                    tmpspec_arr = bkgspec_all_bkgmap[bkgindex0,tmpchmin:tmpchmax]
                    spec_cnt= np.sum(tmpspec_arr)
                    lc_bkg[ii] = lc_bkg[ii] + spec_cnt/np.sum(tmpexpo_arr)
                    #print(lc_time[ii]-lc_time[0],lc_time[ii]-tt0,tt1-lc_time[ii],bkgindex0,lc_bkg[ii])
    
    outname = outnam + '_all.lc'
    write_lcurve(outname,lc_time,lc_bkg,lc_hdr)
    for idid in xrange(0,1):
        print("For detector:", idid)
        for ii in xrange(0,lc_num):
            for jj in xrange(0,bkgmap_num):
                bkgindex0 = idid*bkgmap_num + jj
                tt0 = bkgmap_start_time[jj]
                tt1 = bkgmap_stop_time[jj] + 1
                tmpchmin = chmin + 3
                tmpchmax = chmax + 3
                if (lc_time[ii]>=tt0) & (lc_time[ii]< tt1):
                    tmpexpo_arr = bkgspec_all_bkgmap[bkgindex0,2]
                    tmpspec_arr = bkgspec_all_bkgmap[bkgindex0,tmpchmin:tmpchmax]
                    spec_cnt= np.sum(tmpspec_arr)
                    lc_bkg[ii] = spec_cnt/np.sum(tmpexpo_arr)
                    print(lc_time[ii]-lc_time[0],lc_time[ii]-tt0,tt1-lc_time[ii],bkgindex0,lc_bkg[ii])
        outname = outnam + '_' + str(idid)+'.lc'
        write_lcurve(outname,lc_time,lc_bkg,lc_hdr)
        plt.figure()
        plt.plot(lc_time-lc_time[0],lc_bkg)
        plt.show()
        print(lc_bkg)
'''

print "Finish."





