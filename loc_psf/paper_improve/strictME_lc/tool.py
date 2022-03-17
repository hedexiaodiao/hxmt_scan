import os
import glob
import numpy as np
from astropy.io import fits as pf
import Quat as quat
#from testroll import *
import testroll
from math import fabs,sin,cos,tan,pi,sqrt
from scipy.optimize import curve_fit
from scipy import signal
import matplotlib.pyplot as plt
plt.switch_backend('tkagg')

def Sin(x,a=1.,f=1.,p=1.,c=0.):
    return a*np.sin(2*np.pi/f * x - p)+c

def Sin2(x,a1=1.,f1=1.,p1=1.,c1=0.,a2=1.,f2=1.,p2=1.,c2=0.,a3=1.,c3=1.):
    return Sin(x,a1,f1,p1,c1)*Sin(x,a2,f2,p2,c2)

def CurveFit(func,xsamp,ysamp,X,p0=None,bounds=None):
    if p0==None:
        p0=[1.]*(func.__code__.co_argcount-1)
    if bounds==None:
        bounds = (-np.inf,np.inf)
    popt,pcov=curve_fit(func,xsamp,ysamp,p0=p0,bounds=bounds,sigma=np.sqrt(ysamp+0.75)+1)
    print (popt)
    #print ('COV:',pcov)
    yvals=func(X,*popt)
    return yvals

def FindGap(inptm,gap=1):
    times = np.copy(inptm)
    a=((times[1:]-times[:-1])>gap)
    c=np.nonzero(a)[0]
    potstart = np.r_[0,c+1]
    potstop = np.r_[c,times.shape[0]-1]
    return potstart,potstop 

def peakdel(tm,seqs,interv,limit =100):
    cores = [np.r_[[-1]*interv,[2]*interv,[-1]*interv],np.r_[[1]*interv,[0],[-1]*interv],np.r_[[-1]*interv,[0],[1]*interv]]
    sigs=[]
    gapb,gape = FindGap(tm,5)
    #plt.plot(tm,seqs,'.')
    #print gapb,gape
    for j in cores:
        sig = []
        for i in xrange(len(gapb)):
            ysg = seqs[gapb[i]:gape[i]+1]
            #print ysg,j,np.array(j)
            sig = np.append(sig,signal.convolve(ysg, np.array(j), mode='same'))
            #print 'notiong'
        sigs.append(sig)
        #plt.plot(tm,sig,'.')
    #print 'test'
    jdg = np.where((sigs[0]<10) & np.abs(sigs[1]<10) & np.abs(sigs[2]<10) & (seqs<limit))[0]
    return jdg



def get_rawdata(data_dir, instrument="HE"):
    if instrument == "HE":
        try:filename     = sorted(glob.glob(data_dir + '/HE/*HE-Evt_FFFFFF_V[1-9]*'))[-1]
        except:print("\nERROR: Event file(Evt) not exist...skip this observation\n");
        try:orbitname    = sorted(glob.glob(data_dir + '/ACS/*_Orbit_*V[1-9]*'))[-1]
        except:print("\nERROR: Orbit file(Orbit) not exist...skip this observation\n");
        try:attname      = sorted(glob.glob(data_dir + '/ACS/*_Att_*V[1-9]*'))[-1]
        except:print("\nERROR: Attitude file(Att) not exist...skip this observation\n");
        try:hvfilename   = sorted(glob.glob(data_dir + '/HE/HXMT*HV_FFFFFF*V[1-9]*'))[-1]
        except:print("\nERROR: High Voltage file(HV) not exist...skip this observation\n");
        try:pmfilename   = sorted(glob.glob(data_dir + '/HE/HXMT*PM_FFFFFF*V[1-9]*'))[-1]
        except:print("\nERROR: Particle Monitor file(PM) not exist...skip this observation\n");
        try:deadfilename = sorted(glob.glob(data_dir + '/HE/HXMT*DTime*V[1-9]*'))[-1]
        except:print("\nERROR: Dead time file(DTime) not exist...skip this observation\n");
        try:tempfilename = sorted(glob.glob(data_dir + '/HE/HXMT*TH*V[1-9]*'))[-1]
        except:print("\nERROR: Temperature file(TH) not exist...skip this observation\n");
        try:ehkfilename  = sorted(glob.glob(data_dir + '/AUX/*_EHK_*V[1-9]*'))[-1]
        except:ehkfilename  =''#print("\nERROR: House Keeping file(EHK) not exist...skip this observation\n");
        rawdata_name = ["EVT", "Orbit", "ATT", "HV", "PM", "DTime", "TH", "EHK"]
        rawdata_content = [filename, orbitname, attname, hvfilename, pmfilename, deadfilename, tempfilename, ehkfilename]
        rawdata_dict = dict(zip(rawdata_name, rawdata_content))
        return rawdata_dict
    elif instrument == "ME":
        try:filename     = sorted(glob.glob(data_dir + '/ME/*ME-Evt_FFFFFF_V[1-9]*'))[-1]
        except:print("\nERROR: Event file(Evt) not exist...skip this observation\n");
        try:orbitname    = sorted(glob.glob(data_dir + '/ACS/*_Orbit_*V[1-9]*'))[-1]
        except:print("\nERROR: Orbit file(Orbit) not exist...skip this observation\n");
        try:attname      = sorted(glob.glob(data_dir + '/ACS/*_Att_*V[1-9]*'))[-1]
        except:print("\nERROR: Attitude file(Att) not exist...skip this observation\n");
        try:tempfilename = sorted(glob.glob(data_dir + '/ME/HXMT*TH*V[1-9]*'))[-1]
        except:print("\nERROR: Temperature file(TH) not exist...skip this observation\n");
        try:instatname   = sorted(glob.glob(data_dir + '/ME/HXMT*InsStat*V[1-9]*'))[-1]
        except:print("\nERROR: Instrument Status file(InsStat) not exist...skip this observation\n");
        try:ehkfilename  = sorted(glob.glob(data_dir + '/AUX/*_EHK_*V[1-9]*'))[-1]
        except:print("\nERROR: House Keeping file(EHK) not exist...skip this observation\n");
        rawdata_name = ["EVT", "Orbit", "ATT", "TH", "EHK", "InsStat"]
        rawdata_content = [filename, orbitname, attname, tempfilename, ehkfilename, instatname]
        rawdata_dict = dict(zip(rawdata_name, rawdata_content))
        return rawdata_dict
    elif instrument == "LE":
        try:filename     = sorted(glob.glob(data_dir + '/LE/*LE-Evt_FFFFFF_V[1-9]*'))[-1]
        except:print("\nERROR: Event file(Evt) not exist...skip this observation\n");
        try:orbitname    = sorted(glob.glob(data_dir + '/ACS/*_Orbit_*V[1-9]*'))[-1]
        except:print("\nERROR: Orbit file(Orbit) not exist...skip this observation\n");
        try:attname      = sorted(glob.glob(data_dir + '/ACS/*_Att_*V[1-9]*'))[-1]
        except:print("\nERROR: Attitude file(Att) not exist...skip this observation\n");
        try:tempfilename = sorted(glob.glob(data_dir + '/LE/HXMT*TH*V[1-9]*'))[-1]
        except:print("\nERROR: Temperature file(TH) not exist...skip this observation\n");
        try:instatname   = sorted(glob.glob(data_dir + '/LE/HXMT*InsStat*V[1-9]*'))[-1]
        except:print("\nERROR: Instrument Status file(InsStat) not exist...skip this observation\n");
        try:ehkfilename  = sorted(glob.glob(data_dir + '/AUX/*_EHK_*V[1-9]*'))[-1]
        except:print("\nERROR: House Keeping file(EHK) not exist...skip this observation\n");
        rawdata_name = ["EVT", "Orbit", "ATT", "TH", "EHK", "InsStat"]
        rawdata_content = [filename, orbitname, attname, tempfilename, ehkfilename, instatname]
        rawdata_dict = dict(zip(rawdata_name, rawdata_content))
        return rawdata_dict


def halfbkg(tmary,rtary,erary,qtm,fraglen = 200,step=2):
    xr1 = np.copy(tmary)
    yr1 = np.copy(rtary)
    er1 = np.copy(erary)
    nmb = (yr1.shape[0]-fraglen)/step+1
    lcdic = {}
    
    for idx in range(0,yr1.shape[0]-fraglen+step,step):
        stat = idx
        stop = idx + fraglen if yr1.shape[0] > idx + fraglen else yr1.shape[0]
        #print(idx,stat,stop)
        frag = yr1[stat:stop]
        fragmin = yr1[stat:stop]
        tmmin = xr1[stat:stop]
        ermin = er1[stat:stop]
        polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2,w = 1./ermin))
        #plt.plot(xr1,yr1,'b')
        #plt.plot(xr1[stat:stop],yr1[stat:stop],'r')
        yp=np.in1d(tmmin,qtm)
        fragmin = fragmin[yp]
        tmmin = tmmin[yp]
        ermin = ermin[yp]
        polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2,w = 1./ermin))
        '''
        while 1:
            #plt.plot(tmmin,polmin(tmmin-xr1[0]),'g')
            fragres = -(polmin(tmmin-xr1[0]) - fragmin)
            avg = np.average(fragres**2/ermin**2,weights=1./ermin)*len(fragmin)
            if avg<=len(fragmin)-1-2:break
            fragres=fragres/ermin
            #plt.plot(tmmin[fragres==fragres.max()],fragmin[fragres==fragres.max()],'k.')
            fragmin = fragmin[fragres!=fragres.max()]
            tmmin = tmmin[fragres!=fragres.max()]
            ermin = ermin[fragres!=fragres.max()]
            polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2,w = 1./ermin))
            #plt.plot(tmmin,polmin(tmmin-xr1[0]),'g')
            #plt.pause(0.001)
        #plt.close('all')
        '''
        seg_cut = (stop-stat)//3
        if (stat==0) | (stop>=yr1.shape[0]):
            for ix in range(len(xr1[stat:stop])):
                if str(xr1[stat:stop][ix]) in lcdic:
                    lcdic[str(xr1[stat:stop][ix])].append(polmin(xr1[stat:stop][ix]-xr1[0]))
                else:
                    lcdic[str(xr1[stat:stop][ix])]= [polmin(xr1[stat:stop][ix]-xr1[0])]
        else:
            for ix in range(len(xr1[stat+seg_cut:stop-seg_cut])):
                if str(xr1[stat+seg_cut:stop-seg_cut][ix]) in lcdic:
                    lcdic[str(xr1[stat+seg_cut:stop-seg_cut][ix])].append(polmin(xr1[stat+seg_cut:stop-seg_cut][ix]-xr1[0]))
                else:
                    lcdic[str(xr1[stat+seg_cut:stop-seg_cut][ix])]= [polmin(xr1[stat+seg_cut:stop-seg_cut][ix]-xr1[0])]
    for key,vals in lcdic.items():
        lcdic[key]=np.median(vals)
    tptm = sorted(lcdic.keys())
    tplc = np.array([lcdic[tmid] for tmid in tptm])
    tptms = np.array([float(tmid) for tmid in tptm])
    return tptms,yr1[np.in1d(xr1,tptms)],tplc

def ssbkg(tmary,rtary):
    tm1 = np.copy(tmary)#[20:-20]
    yr1 = np.copy(rtary)#[20:-20]
    xr1 = np.arange(tm1.shape[0])
    xf1 = np.copy(tm1-tm[0])
    interv = np.median(tm[1:]-tm[:-1]);print interv
    #mr1 = signal.medfilt(rtary,31)[20:-20]
    #mr1 = np.copy(rtary)#[20:-20]
    #cores = [np.r_[np.linspace(0,-1,33),[1]*33,np.linspace(-1,0,33)],np.r_[np.linspace(0,1,33),[0],np.linspace(-1,0,49)],np.r_[np.linspace(0,-1,49),[0],np.linspace(1,0,49)]]
    cores = [np.r_[[-1]*int(99/interv),[2]*int(99/interv),[-1]*int(99/interv)],np.r_[[1]*int(99/interv),[0],[-1]*int(99/interv)],np.r_[[-1]*int(99/interv),[0],[1]*int(99/interv)]]
    sigs=[]
    gapb,gape = FindGap(tm1,interv*5)
    for j in cores:
        sig = []
        for i in xrange(len(gapb)):
            ysg = yr1[gapb[i]:gape[i]+1]
            sig = np.append(sig,signal.convolve(ysg, j, mode='same'))
        sigs.append(sig)
    jdg = np.where((sigs[0]<250/interv)&np.abs(sigs[1]<199/interv)&np.abs(sigs[2]<199/interv)&(yr1<55*interv))[0]
    jdgx = xr1[jdg]
    jdgy  = yr1[jdg]
    plt.plot(tm1,yr1)
    plt.plot(tm1[jdgx],jdgy,'.')
    gapb,gape = FindGap(tm1,1000/interv)
    print gapb,gape
    for i in xrange(len(gapb)):
        xb1 = np.intersect1d(jdgx,xr1[gapb[i]:gape[i]+1])
        yb1 = yr1[xb1]
        print xb1[0],xb1[-1]
        z=CurveFit(Sin2,(xf1[xb1]-xf1[xb1][0])/interv,yb1,(xf1[xb1]-xf1[xb1][0])/interv,p0=[1.,3001/interv,1,0,1,2001/interv,1,0,1,0],\
                   bounds=([-np.inf,1000/interv,-np.inf,-np.inf,-np.inf,1000/interv,-np.inf,-np.inf,-np.inf,-np.inf],[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]))
        z=CurveFit(Sin2,(xf1[xb1]-xf1[xb1][0])[z<30*interv]/interv,yb1[z<30*interv],(xf1[xb1[0]:xb1[-1]]-xf1[xb1[0]:xb1[-1]][0])/interv,p0=[1.,3001/interv,1,0,1,2001/interv,1,0,1,0],\
                   bounds=([-np.inf,1000/interv,-np.inf,-np.inf,-np.inf,1000/interv,-np.inf,-np.inf,-np.inf,-np.inf],[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]))
        x_r,y_r,bkg_r = tm1[xb1[0]:xb1[-1]],yr1[xb1[0]:xb1[-1]],z
        if i ==0:
            bkg = bkg_r
            xr_sel = x_r
            yr_sel = y_r
        if i >0:
            try:
                bkg = np.r_[bkg,bkg_r]
                xr_sel = np.r_[xr_sel,x_r]
                yr_sel = np.r_[yr_sel,y_r]
            except UnboundLocalError:
                bkg = bkg_r
                xr_sel = x_r
                yr_sel = y_r
        plt.plot(x_r,bkg_r)
    plt.show()
    return xr_sel,yr_sel,bkg

def bkg_new_ary(x1,y1,err1,qtm,fraglen = 200,step=2):
    xb = x1.copy()
    yb = y1.copy()
    errb=err1.copy()
    Ygti = ((xb[1:]-xb[:-1])>100)
    c=np.nonzero(Ygti)[0]
    pot=np.r_[0,c+1,Ygti.shape[0]]
    print pot
    for i in xrange(len(pot)-1): 
        try:
            yb1 = yb[pot[i]:pot[i+1]][2:-2]
            xb1 = xb[pot[i]:pot[i+1]][2:-2]
            erb1= errb[pot[i]:pot[i+1]][2:-2]
        except IndexError:
            raw_input("Indexerror ,recompose you code!")
            yb1 = yb[pot[i]:pot[i+1]]
            xb1 = xb[pot[i]:pot[i+1]]
        if yb1.shape[0]>=400:
            x_r,y_r,bkg_r = halfbkg(xb1,yb1,erb1,qtm,fraglen,step)
            print len(x_r),len(y_r),len(bkg_r)
        else:
            continue
        if i ==0:
            bkg = bkg_r
            xr_sel = x_r
            yr_sel = y_r
        if i >0:
            try:
                bkg = np.r_[bkg,bkg_r]
                xr_sel = np.r_[xr_sel,x_r]
                yr_sel = np.r_[yr_sel,y_r]
            except UnboundLocalError:
                bkg = bkg_r
                xr_sel = x_r
                yr_sel = y_r
    return xr_sel,yr_sel,bkg


def bkg_new(lcfile):
    hdu = pf.open(lcfile)
    tb = hdu[1].data
    x = tb.field(0)
    y = tb.field(1)
    err=tb.field(2)
    xr_sel,yr_sel,bkg=bkg_new_ary(x,y,err)
    hdu.close()
    del hdu
    return xr_sel,yr_sel,bkg


def psfpeak(attfile,ra,dec,instr='ME',box=0):
    if instr=='ME':flag=1
    if instr=='HE':flag=0
    if instr=='LE':flag=2
    atthd = pf.open(attfile)
    atdt = atthd[3].data
    qtime = atdt.field(0)
    q1_list = atdt.field(1)
    q2_list = atdt.field(2)
    q3_list = atdt.field(3)
    atthd.close();del atthd
    Talfa_bound = np.array([5.7,4.0,6.0])
    Tbeta_bound = np.array([1.1,1.0,1.6])
    dcm_b_f = [np.array([[0,1,0],[0,0,1],[1,0,0]]),np.array([[0,0,-1],[0,1,0],[1,0,0]]),np.array([[0,0,1],[0,1,0],[-1,0,0]])]
    mtx_list=[]
    for i in range(0,q1_list.shape[0]):
        quat1 = [q1_list[i],q2_list[i],q3_list[i]]
        mtx_list.append(testroll.quat2mtx(quat1))
    
    D_alfa = []
    D_beta= []
    paras0 = np.loadtxt('/sharefs/hbkg/user/nangyi/lick_tools/para_files/Parabolic_%s_%s_201904.txt'%('HE','all'))
    paras1 = np.loadtxt('/sharefs/hbkg/user/nangyi/lick_tools/para_files/Parabolic_%s_%s_201904.txt'%('ME','all'))
    paras2 = np.loadtxt('/sharefs/hbkg/user/nangyi/lick_tools/para_files/Parabolic_%s_%s_201904.txt'%('LE','all'))
    Tpsai=[paras0[::9,0],paras1[::9,0],paras2[::9,0]]
    Ttheta=[paras0[1::9,0],paras1[1::9,0],paras2[1::9,0]]
    Tphi=[paras0[2::9,0],paras1[2::9,0],paras2[2::9,0]]
    abound = Talfa_bound + np.array([paras0[3::9,0],paras1[3::9,0],paras2[3::9,0]])
    bbound = Tbeta_bound + np.array([paras0[4::9,0],paras1[4::9,0],paras2[4::9,0]])
    Tpa,Tpb,Tpc,Tpd = [paras0[5::9,0],paras1[5::9,0],paras2[5::9,0]],[paras0[6::9,0],paras1[6::9,0],paras2[6::9,0]],[paras0[7::9,0],paras1[7::9,0],paras2[7::9,0]],[paras0[8::9,0],paras1[8::9,0],paras2[8::9,0]]
    flux=np.zeros(len(mtx_list))
    PI=3.14159265358979323846
    roll=box/180.0*PI
    roll_indx = int(-box/60+1)
    psai=Tpsai[flag][roll_indx]
    theta=Ttheta[flag][roll_indx]
    phi=Tphi[flag][roll_indx]
    alfa_bound = abound[flag][roll_indx]
    beta_bound = bbound[flag][roll_indx]
    pa,pb,pc,pd =Tpa[flag][roll_indx],Tpb[flag][roll_indx],Tpc[flag][roll_indx],Tpd[flag][roll_indx]
    s_step=0
    Nflux = len(mtx_list)
    rot_quat =  quat.Quat((psai,theta,phi))
    x,y,z = testroll.sph2cart(ra*np.pi/180,dec*np.pi/180,1)
    xyz = [x,y,z]
    dcm_b_f_rot = rot_quat.transform.T
    for s_step in range(0,Nflux):
        mtx = mtx_list[s_step]
        delta_alfa0,delta_beta0,xr,yr,zr = testroll.quat2delta2_rot(mtx,dcm_b_f[flag],xyz,dcm_b_f_rot)
        cosz = np.sqrt(zr**2/(xr**2+yr**2+zr**2))
        delta_alfa=fabs(delta_alfa0*cos(roll)-delta_beta0*sin(roll))
        delta_beta=fabs(delta_alfa0*sin(roll)+delta_beta0*cos(roll))
        if delta_alfa<alfa_bound and delta_beta<beta_bound:
            factor=(1.0-tan(delta_alfa/180.0*PI)/tan(alfa_bound/180.0*PI))*\
            (1.0-tan(delta_beta/180.0*PI)/tan(beta_bound/180.0*PI))
        else:
            factor=0
        factor= factor/(sqrt(tan(delta_alfa/180.0*PI)*tan(delta_alfa/180.0*PI)+tan(delta_beta/180.0*PI)*tan(delta_beta/180*PI)+1))
        rt = pa * delta_alfa**2 * delta_beta**2 + pb*delta_alfa**2 +pc*delta_beta**2 + pd
        flux[s_step]=factor*rt
    return flux

class Psfmodel(object):
    def __init__(self, attfile):
        atthd = pf.open(attfile)
        atdt = atthd[3].data
        qtime = atdt.field(0)
        q1_list = atdt.field(1)
        q2_list = atdt.field(2)
        q3_list = atdt.field(3)
        atthd.close();del atthd
        Talfa_bound = np.array([5.7,4.0,6.0])
        Tbeta_bound = np.array([1.1,1.0,1.6])
        dcm_b_f = [np.array([[0,1,0],[0,0,1],[1,0,0]]),np.array([[0,0,-1],[0,1,0],[1,0,0]]),np.array([[0,0,1],[0,1,0],[-1,0,0]])]
        mtx_list=[]
        for i in range(0,q1_list.shape[0]):
            quat1 = [q1_list[i],q2_list[i],q3_list[i]]
            mtx_list.append(testroll.quat2mtx(quat1))
        self.mtxs  = mtx_list
        self.qtime = qtime
    def Flux(self,ra,dec,instr='ME',box=0):
        if instr=='ME':flag=1
        if instr=='HE':flag=0
        if instr=='LE':flag=2
        flux=np.zeros(len(self.mtxs))
        Talfa_bound = np.array([5.7,4.0,6.0])
        Tbeta_bound = np.array([1.1,1.0,1.6])
        dcm_b_f = [np.array([[0,1,0],[0,0,1],[1,0,0]]),np.array([[0,0,-1],[0,1,0],[1,0,0]]),np.array([[0,0,1],[0,1,0],[-1,0,0]])]
        paras0 = np.loadtxt('/sharefs/hbkg/user/nangyi/lick_tools/para_files/Parabolic_%s_%s_201904.txt'%('HE','all'))
        paras1 = np.loadtxt('/sharefs/hbkg/user/nangyi/lick_tools/para_files/Parabolic_%s_%s_201904.txt'%('ME','all'))
        paras2 = np.loadtxt('/sharefs/hbkg/user/nangyi/lick_tools/para_files/Parabolic_%s_%s_201904.txt'%('LE','all'))
        Tpsai=[paras0[::9,0],paras1[::9,0],paras2[::9,0]]
        Ttheta=[paras0[1::9,0],paras1[1::9,0],paras2[1::9,0]]
        Tphi=[paras0[2::9,0],paras1[2::9,0],paras2[2::9,0]]
        abound = Talfa_bound + np.array([paras0[3::9,0],paras1[3::9,0],paras2[3::9,0]])
        bbound = Tbeta_bound + np.array([paras0[4::9,0],paras1[4::9,0],paras2[4::9,0]])
        Tpa,Tpb,Tpc,Tpd = [paras0[5::9,0],paras1[5::9,0],paras2[5::9,0]],[paras0[6::9,0],paras1[6::9,0],paras2[6::9,0]],[paras0[7::9,0],paras1[7::9,0],paras2[7::9,0]],[paras0[8::9,0],paras1[8::9,0],paras2[8::9,0]]
        roll=box/180.0*pi
        roll_indx = int(-box/60+1)
        psai=Tpsai[flag][roll_indx]
        theta=Ttheta[flag][roll_indx]
        phi=Tphi[flag][roll_indx]
        alfa_bound = abound[flag][roll_indx]
        beta_bound = bbound[flag][roll_indx]
        pa,pb,pc,pd =Tpa[flag][roll_indx],Tpb[flag][roll_indx],Tpc[flag][roll_indx],Tpd[flag][roll_indx]
        rot_quat =  quat.Quat((psai,theta,phi))
        x,y,z = testroll.sph2cart(ra*pi/180,dec*pi/180,1)
        xyz = [x,y,z]
        dcm_b_f_rot = rot_quat.transform
        for s_step in range(0,len(self.mtxs)):
            mtx = self.mtxs[s_step]
            delta_alfa0,delta_beta0,xr,yr,zr = testroll.quat2delta2_rot(mtx,dcm_b_f[flag],xyz,dcm_b_f_rot)
            cosz = np.sqrt(zr**2/(xr**2+yr**2+zr**2))
            delta_alfa=fabs(delta_alfa0*cos(roll)-delta_beta0*sin(roll))
            delta_beta=fabs(delta_alfa0*sin(roll)+delta_beta0*cos(roll))
            if delta_alfa<alfa_bound and delta_beta<beta_bound:
                factor=(1.0-tan(delta_alfa/180.0*pi)/tan(alfa_bound/180.0*pi))*\
                (1.0-tan(delta_beta/180.0*pi)/tan(beta_bound/180.0*pi))
            else:
                factor=0
            factor= factor/(sqrt(tan(delta_alfa/180.0*pi)*tan(delta_alfa/180.0*pi)+tan(delta_beta/180.0*pi)*tan(delta_beta/180*pi)+1))
            rt = pa * delta_alfa**2 * delta_beta**2 + pb*delta_alfa**2 +pc*delta_beta**2 + pd
            flux[s_step]=factor*rt
        return flux


def Gauss(x,*param):
    return param[0]*np.exp(-np.power(x - param[1], 2.) / (2 * np.power(param[2], 2.)))

def Gaussfit(x,y,p0=[1,1,1]):
    popt,pcov = curve_fit(Gauss,x,y,p0=p0)
    sigma = np.sqrt(np.diag(pcov))
    y1=Gauss(x,*popt)
    return y1,popt,sigma


def genmebd(screenfile,outfile,bdfile=None):
    hd=pf.open(screenfile)
    det = hd[1].data['Det_ID']
    det=det[det<1728]
    a,b = np.histogram(det,bins=np.arange(-0.5,576,1.))
    did = np.arange(576)
    y,z = np.histogram(a,bins=100)
    x=(z[1:]+z[:-1])/2.
    c,p,s=Gaussfit(x,y,p0=[10,np.median(a),np.sqrt(np.median(a))])
    plt.plot(x,y);plt.plot(x,c)
    plt.show()
    badid = did[(a<p[1]-2*p[2])|(a>p[1]+2*p[2])]
    a,b = np.histogram(det,bins=np.arange(576-0.5,1152,1.))
    did = np.arange(576,1152)
    y,z = np.histogram(a,bins=100)
    x=(z[1:]+z[:-1])/2.
    c,p,s=Gaussfit(x,y,p0=[10,np.median(a),np.sqrt(np.median(a))])
    plt.plot(x,y);plt.plot(x,c)
    plt.show()
    badid = np.append(did[(a<p[1]-2*p[2])|(a>p[1]+2*p[2])],badid)
    a,b = np.histogram(det,bins=np.arange(1152-0.5,1728,1.))
    did = np.arange(1152,1728)
    y,z = np.histogram(a,bins=100)
    x=(z[1:]+z[:-1])/2.
    c,p,s=Gaussfit(x,y,p0=[10,np.median(a),np.sqrt(np.median(a))])
    plt.plot(x,y);plt.plot(x,c)
    plt.show()
    badid = np.append(did[(a<p[1]-2*p[2])|(a>p[1]+2*p[2])],badid)
    hd.close();del hd
    badid=set(badid)
    if bdfile is not None:
        hd=pf.open(bdfile)
        oldbd=hd[1].data.field(0)
        badid = badid.union(oldbd)
        hd.close();del hd
    badnum=list(badid)
    hd=pf.open('/sharefs/hbkg/data/SCAN/nangyi/ME_data/gti/mebadall.fits')
    prihdr=hd[1].header
    lth = len(badnum)
    col1 = pf.Column(name = 'DetID', format='I',array = np.sort(badnum))
    col2 = pf.Column(name = 'TIMERANGE', format='20A',array = np.array(['0']*lth))
    col3 = pf.Column(name = 'TIMERANGE2', format='20A',array = np.array(['INDEF']*lth))
    col4 = pf.Column(name = 'TYPE', format='20A',array = np.array(['Bad']*lth))
    col5 = pf.Column(name = 'STATUS', format='B',array = np.array([0]*lth))
    cols = pf.ColDefs([col1,col2,col3,col4,col5])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    pf.conf.extension_name_case_sensitive=True
    tbhdu.name='detectorStatus'
    tbhdu.header = prihdr
    tbhdu.writeto(outfile,overwrite=True)
    return list(badid)

