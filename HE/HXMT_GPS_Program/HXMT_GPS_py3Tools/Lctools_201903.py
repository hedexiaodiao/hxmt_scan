#!/hxmt/soft/Develop/anaconda2/bin/python
import numpy as np
from astropy.io import fits as pf
import sys,os
from glob import glob as glob
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import Findfile
from Astrotools import Sel_Median
run = Findfile.findfile()

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

def FindGap(inptm,gap=1):
    times = np.copy(inptm)
    a=((times[1:]-times[:-1])>gap)
    c=np.nonzero(a)[0]
    potstart = np.r_[0,c+1]
    potstop = np.r_[c,times.shape[0]-1]
    return potstart,potstop

def GenArraFromGap(tst,tsp,gap=1):
    tmlt = np.array([])
    for ix,it in enumerate(tst):
        try:
            tmlt=np.append(tmlt,np.arange(it,tsp[ix]+1,gap))
        except:
            pass
    return tmlt

def Segments_Gap(inptm,gap=1):
    index = range(len(inptm))
    potstart,potstop = FindGap(inptm,gap=gap)
    segs = []
    ind = []
    for i in range(len(potstart)):
        if potstart[i]<potstop[i]:
            segs.append((inptm[potstart[i] : potstop[i]+1]).tolist())
            ind.append(index[potstart[i] : potstop[i]+1])
    return segs,ind

def FindLC_basefile(key,fpath,INST='LE',EHKFILE=True):
    '''
    return order:
    EVTfile,THfile,ATTfile,(LE:INSfile//HE:DTfile,HVfile,PMfile),EHKfile//ORTfile
    '''
    shortkey = key[:-2]
    returnlist = []
    evtlisttp1 = sorted(run(key+"_%s-Evt_FFFFFF_V"%INST,fpath))
    print(fpath)
    print(evtlisttp1)
    try:
        float(evtlisttp1[-1].split('/')[-1][-10])
        evtfile = evtlisttp1[-1]
    except:
        evtfile = evtlisttp1[-2]
    returnlist.append(evtfile)
    tplisttp1 = sorted(run(key+"_%s-TH_FFFFFF_V"%INST,fpath))
    try:
        float(tplisttp1[-1].split('/')[-1][-10])
        tpfile = tplisttp1[-1]
    except:
        tpfile = tplisttp1[-2]
    returnlist.append(tpfile)
    attlisttp1 = sorted(run(key+"_Att_FFFFFF_V",fpath))
    if len(attlisttp1)==0:
        attlisttp1 = sorted(run(shortkey+"_Att_FFFFFF_V",fpath))
    try:
        float(attlisttp1[-1].split('/')[-1][-10])
        attfile = attlisttp1[-1]
    except:
        attfile = attlisttp1[-2]
    returnlist.append(attfile)
    if INST=='LE':
        Inslisttp1 = sorted(run(key+"_%s-InsStat_FFFFFF_V"%INST,fpath))
        try:
            float(Inslisttp1[-1].split('/')[-1][-10])
            Insfile = Inslisttp1[-1]
        except:
            Insfile = Inslisttp1[-2]
        returnlist.append(Insfile)
    elif INST=='ME':
        pass
    elif INST=='HE':
        dtlisttp1 = sorted(run(key+"_HE-DTime_FFFFFF_V",fpath))
        try:
            float(dtlisttp1[-1].split('/')[-1][-10])
            dtfile = dtlisttp1[-1]
        except:
            dtfile = dtlisttp1[-2]
        returnlist.append(dtfile)
        hvlisttp1 = sorted(run(key+"_HE-HV_FFFFFF_V",fpath))
        try:
            float(hvlisttp1[-1].split('/')[-1][-10])
            hvfile = hvlisttp1[-1]
        except:
            hvfile = hvlisttp1[-2]
        returnlist.append(hvfile)
        pmlisttp1 = sorted(run(key+"_HE-PM_FFFFFF_V",fpath))
        try:
            float(pmlisttp1[-1].split('/')[-1][-10])
            pmfile = pmlisttp1[-1]
        except:
            pmfile = pmlisttp1[-2]
        returnlist.append(pmfile)
    if EHKFILE:
        ehklisttp1 = sorted(run(key+"_EHK_FFFFFF_V",fpath))
        if len(ehklisttp1)==0:
            ehklisttp1 = sorted(run(shortkey+"_EHK_FFFFFF_V",fpath))
        try:
            float(ehklisttp1[-1].split('/')[-1][-10])
            ehkfile = ehklisttp1[-1]
        except:
            ehkfile = ehklisttp1[-2]
        returnlist.append(ehkfile)
    else:
        ortlisttp1 = sorted(run(key +"_Orbit_FFFFFF_V",fpath))
        if len(ortlisttp1)==0:
            ortlisttp1 = sorted(run(shortkey +"_Orbit_FFFFFF_V",fpath))
        try:
            float(ortlisttp1[-1].split('/')[-1][-10])
            ortfile = ortlisttp1[-1]
        except:
            ortfile = ortlisttp1[-2]
        returnlist.append(ortfile)
    return returnlist

def he_nozero_pot(pifile,ehkfile):
    hd  = pf.open(pifile)
    ehk = pf.open(ehkfile)

    evtdat = hd[1].data
    timm = evtdat.field(0)
    chnl = evtdat.field(2)
    pulw = evtdat.field(3)
    acdd = evtdat.field(4)
    etyp = evtdat.field(5)
    d_id = evtdat.field(1)
    pi = evtdat.field(7)

    ehkdat = ehk[1].data
    ehktm= ehkdat.field(0)
    #elv = ehkdat.field(10)
    #cor = ehkdat.field('COR')
    #lat = ehkdat.field(9)
    #delv=ehkdat.field("DYE_ELV")
    tsaa=ehkdat.field("T_SAA")
    nsaa=ehkdat.field("TN_SAA")

    acdyn= np.sum((acdd == 0),axis=1)
    Yseq =(pi <= 87)&(pi >0)&(etyp == 0)&(acdyn == 18)&(pulw < 71)&(pulw > 29)
    timm_sel = timm[(Yseq)]
    d_id_sel = d_id[(Yseq)]

    tstart = np.ceil(min(timm_sel)) #+ ehktm[0] - np.floor(ehktm[0])

    tstop  = np.floor(max(timm_sel)) #+ ehktm[0] - np.floor(ehktm[0])
    print(tstart)
    print(tstop)
    tm_id1 = timm_sel[( (d_id_sel == 1 )|(d_id_sel == 3)|(d_id_sel == 8)|(d_id_sel == 11)|(d_id_sel == 12) )]
    tm_id0 = timm_sel[( (d_id_sel == 0 )|(d_id_sel == 4)|(d_id_sel == 7)|(d_id_sel == 14)|(d_id_sel == 17) )]
    tm_id5 = timm_sel[( (d_id_sel == 5 )|(d_id_sel == 6)|(d_id_sel == 10)|(d_id_sel == 13)|(d_id_sel == 15) )]
    #tm_id16 = timm_sel[(d_id_sel == 16)]
    #tm_id2 = timm_sel[(d_id_sel == 2)]
    #tm_id9 = timm_sel[(d_id_sel == 9)]
    tm_binset_1d = np.arange(tstart+0.5,tstop+1,1)
    phon_id0,tbins_id0 = np.histogram(tm_id0,tm_binset_1d)
    phon_id1,tbins_id1 = np.histogram(tm_id1,tm_binset_1d)
    phon_id5,tbins_id5 = np.histogram(tm_id5,tm_binset_1d)
    #phon_id16,tbins_id16 = np.histogram(tm_id16,tm_binset_1d)
    #phon_id2,tbins_id2 = np.histogram(tm_id2,tm_binset_1d)
    #phon_id9,tbins_id9 = np.histogram(tm_id9,tm_binset_1d)

    times=ehktm[((tsaa>0)&(nsaa>0))]
    timeid0 = tbins_id0[:-1][(phon_id0>0)]
    timeid1 = tbins_id1[:-1][(phon_id1>0)]
    timeid5 = tbins_id5[:-1][(phon_id5>0)]
    Yseq = (etyp == 1)&(acdyn == 18)
    timm_sel = timm[(Yseq)]
    d_id_sel = d_id[(Yseq)]
    csi1d,bins = np.histogram(timm_sel,tm_binset_1d)
    print('normal time:::',(csi1d>0).sum())
    if (csi1d>0).sum()>600:
        flag=1
        sm_rad = smooth(csi1d,51)
        xcenter = bins[:-1][sm_rad>0]
        med_sel = Sel_Median(-sm_rad[sm_rad>0],3.5)
        tm_medsel = xcenter[med_sel[1]]
        tmsl = np.array(Uniform([timeid0.tolist(),timeid1.tolist(),timeid5.tolist(),times.tolist(),tm_medsel.tolist()]))
    else:
        flag=0
        sm_rad = smooth(csi1d,51)
        xcenter = bins[:-1].copy()
        med_sel = (sm_rad==0)
        tm_medsel = xcenter[med_sel]
        tmsl = np.array(Uniform([timeid0.tolist(),timeid1.tolist(),timeid5.tolist(),times.tolist(),tm_medsel.tolist()]))
    tmsl =(tmsl+0.5).astype(np.int)
    b = ((tmsl[1:]-tmsl[:-1])>1)
    c = np.nonzero(b)[0]
    potstart = np.r_[0,c+1]
    potstop = np.r_[c,tmsl.shape[0]-1]
    timestart= tmsl[potstart]
    timestop = tmsl[potstop]
    fig = plt.figure(figsize=(12,10))
    plt.plot((bins[1:]+bins[:-1])/2.,csi1d)
    plt.title('Evtype == 1',fontsize=10)
    plt.xlabel('Time /s',fontsize=10)
    plt.autoscale(enable=True, tight=True)
    figname = pifile.split('pi.f')[0]
    plt.savefig(figname + 'radiation.png')
    plt.close('all')
    hd.close()
    ehk.close()
    del hd
    del ehk
    return [timestart,timestop],flag

def he_nozero_pot_CSI(pifile,ehkfile):
    hd  = pf.open(pifile)
    ehk = pf.open(ehkfile)

    evtdat = hd[1].data
    timm = evtdat.field(0)
    chnl = evtdat.field(2)
    pulw = evtdat.field(3)
    acdd = evtdat.field(4)
    etyp = evtdat.field(5)
    d_id = evtdat.field(1)
    pi = evtdat.field(7)

    ehkdat = ehk[1].data
    ehktm= ehkdat.field(0)
    #elv = ehkdat.field(10)
    #cor = ehkdat.field('COR')
    #lat = ehkdat.field(9)
    #delv=ehkdat.field("DYE_ELV")
    tsaa=ehkdat.field("T_SAA")
    nsaa=ehkdat.field("TN_SAA")

    acdyn= np.sum((acdd == 0),axis=1)
    Yseq =(pi <= 87)&(pi >0)&(etyp == 0)&(acdyn == 18)&(pulw < 71)&(pulw > 29)
    timm_sel = timm[(Yseq)]
    d_id_sel = d_id[(Yseq)]

    tstart = np.ceil(min(timm_sel)) #+ ehktm[0] - np.floor(ehktm[0])

    tstop  = np.floor(max(timm_sel)) #+ ehktm[0] - np.floor(ehktm[0])
    print(tstart)
    print(tstop)
    tm_id1 = timm_sel[( (d_id_sel == 1 )|(d_id_sel == 3)|(d_id_sel == 8)|(d_id_sel == 11)|(d_id_sel == 12) )]
    tm_id0 = timm_sel[( (d_id_sel == 0 )|(d_id_sel == 4)|(d_id_sel == 7)|(d_id_sel == 14)|(d_id_sel == 17) )]
    tm_id5 = timm_sel[( (d_id_sel == 5 )|(d_id_sel == 6)|(d_id_sel == 10)|(d_id_sel == 13)|(d_id_sel == 15) )]
    #tm_id16 = timm_sel[(d_id_sel == 16)]
    #tm_id2 = timm_sel[(d_id_sel == 2)]
    #tm_id9 = timm_sel[(d_id_sel == 9)]
    tm_binset_1d = np.arange(tstart+0.5,tstop+1,1)
    phon_id0,tbins_id0 = np.histogram(tm_id0,tm_binset_1d)
    phon_id1,tbins_id1 = np.histogram(tm_id1,tm_binset_1d)
    phon_id5,tbins_id5 = np.histogram(tm_id5,tm_binset_1d)
    #phon_id16,tbins_id16 = np.histogram(tm_id16,tm_binset_1d)
    #phon_id2,tbins_id2 = np.histogram(tm_id2,tm_binset_1d)
    #phon_id9,tbins_id9 = np.histogram(tm_id9,tm_binset_1d)

    times=ehktm[((tsaa>0)&(nsaa>0))]
    timeid0 = tbins_id0[:-1][(phon_id0>0)]+0.5
    timeid1 = tbins_id1[:-1][(phon_id1>0)]+0.5
    timeid5 = tbins_id5[:-1][(phon_id5>0)]+0.5
    print(timeid0)
    Yseq = (etyp == 1)&(acdyn == 18)
    timm_sel = timm[(Yseq)]
    d_id_sel = d_id[(Yseq)]
    csi1d,bins = np.histogram(timm_sel,tm_binset_1d)
    xcenter = (bins[1:]+bins[:-1])/2
    print(xcenter, xcenter.shape)
    print('normal time:::',(csi1d>0).sum())
    if (csi1d>0).sum()>600:
        flag=1
        sm_rad = np.zeros(len(csi1d))
        tst,tsp =FindGap(xcenter[csi1d>0])
        print('tst,tsp',tst,tsp)
        for ix,it in enumerate(tst):
            start = np.arange(len(csi1d))[csi1d>0][it]
            stop = np.arange(len(csi1d))[csi1d>0][tsp[ix]]
            print('st',start,start.shape,stop,stop.shape)
            try:
                sm_seg = smooth(csi1d[start:stop+1],51)
                print('try',sm_seg,sm_seg.shape)
            except:
                sm_seg = np.zeros(stop-start+1)
                print('except', sm_seg, sm_seg.shape)
            sm_rad[start:stop+1] = sm_seg[:]
        print('sm_rad',sm_rad,sm_rad.shape,sm_rad>300,(sm_rad>300).shape,np.sum((np.array(sm_rad>300)).astype(int)))
        med_sel = sm_rad>300
        print(med_sel,med_sel[1],med_sel.shape)
        tm_medsel = xcenter[med_sel]
        print('if',timeid0.shape,timeid1.shape,timeid5.shape,times.shape,tm_medsel.shape)
        tmsl = np.array(Uniform([timeid0.tolist(),timeid1.tolist(),timeid5.tolist(),times.tolist(),tm_medsel.tolist()]))
    else:
        flag=0
        sm_rad = smooth(csi1d,51)
        med_sel = (sm_rad==0)
        tm_medsel = xcenter[med_sel]
        print('else', timeid0.shape, timeid1.shape, timeid5.shape, times.shape, tm_medsel.shape)
        tmsl = np.array(Uniform([timeid0.tolist(),timeid1.tolist(),timeid5.tolist(),times.tolist(),tm_medsel.tolist()]))
    tmsl =tmsl.astype(np.int)
    b = ((tmsl[1:]-tmsl[:-1])>1)
    c = np.nonzero(b)[0]
    potstart = np.r_[0,c+1]
    potstop = np.r_[c,tmsl.shape[0]-1]
    print(potstart,potstop)
    print(tmsl.shape)
    timestart= tmsl[potstart]
    timestop = tmsl[potstop]
    print(timestart)
    print(timestop)
    sel_start=[]
    sel_stop=[]
    for i in range(len(timestart)):
        if (timestop[i]-timestart[i])>10:
            sel_start.append(timestart[i])
            sel_stop.append(timestop[i])
    fig = plt.figure(figsize=(12,10))
    plt.plot(xcenter,csi1d)
    print(xcenter[0],tmsl[0])
    plt.plot(tmsl,csi1d[np.in1d(xcenter,tmsl)],'.')
    plt.title('Evtype == 1',fontsize=10)
    plt.xlabel('Time /s',fontsize=10)
    plt.autoscale(enable=True, tight=True)
    figname = pifile.split('pi.f')[0]
    plt.savefig(figname + 'radiation.png')
    plt.close('all')
    hd.close()
    ehk.close()
    del hd
    del ehk
    return [sel_start,sel_stop],flag

def norm_18lc_v2(key,lcpath):
    tm0 = pf.open(lcpath + key + "_g0_0.lc")[1].data.field(0)
    tm1 = pf.open(lcpath + key + "_g0_1.lc")[1].data.field(0)
    tm2 = pf.open(lcpath + key + "_g0_2.lc")[1].data.field(0)
    tm3 = pf.open(lcpath + key + "_g0_3.lc")[1].data.field(0)
    tm4 = pf.open(lcpath + key + "_g0_4.lc")[1].data.field(0)
    tm5 = pf.open(lcpath + key + "_g0_5.lc")[1].data.field(0)
    tm6 = pf.open(lcpath + key + "_g0_6.lc")[1].data.field(0)
    tm7 = pf.open(lcpath + key + "_g0_7.lc")[1].data.field(0)
    tm8 = pf.open(lcpath + key + "_g0_8.lc")[1].data.field(0)
    tm9 = pf.open(lcpath + key + "_g0_9.lc")[1].data.field(0)
    tm10 = pf.open(lcpath + key + "_g0_10.lc")[1].data.field(0)
    tm11 = pf.open(lcpath + key + "_g0_11.lc")[1].data.field(0)
    tm12 = pf.open(lcpath + key + "_g0_12.lc")[1].data.field(0)
    tm13 = pf.open(lcpath + key + "_g0_13.lc")[1].data.field(0)
    tm14 = pf.open(lcpath + key + "_g0_14.lc")[1].data.field(0)
    tm15 = pf.open(lcpath + key + "_g0_15.lc")[1].data.field(0)
    tm16 = pf.open(lcpath + key + "_g0_16.lc")[1].data.field(0)
    tm17 = pf.open(lcpath + key + "_g0_17.lc")[1].data.field(0)

    tmsl = Uniform([tm0.tolist(),tm1.tolist(),tm2.tolist(),tm3.tolist(),tm4.tolist(),tm5.tolist(),tm6.tolist(),tm7.tolist(),tm8.tolist(),tm9.tolist(),tm10.tolist(),tm11.tolist(),tm12.tolist(),tm13.tolist(),tm14.tolist(),tm15.tolist(),tm16.tolist(),tm17.tolist()])
    tmsl = np.array(tmsl)
    for i in range(18):
        filename = lcpath + key + "_g0_%s.lc"%i
        hd = pf.open(filename,mode='update')
        tb = hd[1].data
        hd[1].data = tb[np.in1d(tb.field(0),tmsl)]
        hd.flush()
        hd.close()
        del hd,tb
    print("Normal Lc OK!")

def norm_3netlc(lclist,lcpath):
    tm0 = pf.open(lclist[0])[1].data.field(0)
    tm1 = pf.open(lclist[1])[1].data.field(0)
    tm2 = pf.open(lclist[2])[1].data.field(0)
    tmsl = Uniform([tm0.tolist(),tm1.tolist(),tm2.tolist()])
    for i in range(3):
        hd = pf.open(lclist[i],mode='update')
        tb = hd[1].data 
        hd[1].data = tb[np.in1d(tb.field(0),tmsl)]
        hd.flush()
        hd.close()
        del hd,tb
    print("Normal Lc OK!")

def norm_3lc_bkg_v2(key,bkgpath):
    tm0 = pf.open(bkgpath + key +'he_netlc_b0.fits')[1].data.field(0)###
    tm1 = pf.open(bkgpath + key +'he_netlc_b1.fits')[1].data.field(0)
    tm2 = pf.open(bkgpath + key +'he_netlc_b2.fits')[1].data.field(0)
    tmsl = Uniform([tm0.tolist(),tm1.tolist(),tm2.tolist()])
    for i in range(3):
        if i>1:
            i=i###
        filename = glob(bkgpath+key+'he_netlc_b%s.fits'%(i))[0]# bkgpath + key + "lc_%s.lc"%i
        print(filename)
        hd = pf.open(filename,mode='update')
        tb = hd[1].data 
        hd[1].data = tb[np.in1d(tb.field(0),tmsl)]
        hd.flush()
        hd.close()
        del hd,tb
    print("Normal Bkg-Fits OK!")

def togeth_18lc_v2(key,lcpath):
    norm_18lc_v2(key,lcpath)
    print("Together LC")
    infilekey= key
    lcid0=[0,4,7,14,17]
    lcid1=[1,3,8,11,12]
    lcid2=[5,6,10,13,15]
    lc0=np.array([])
    lc2=np.array([])
    lc1=np.array([])
    for i in range(18):
        hd=pf.open(lcpath+infilekey+"_g0_%s.lc"%i)
        tb = hd[1].data
        times = tb.field(0)
        allcts = tb.field(1)
        if (i!=0)&(i!=1)&(i!=5):
            if i in lcid0:
                lc0=lc0 + allcts
                #print allcts[0],allcts[1],allcts[2],"lc0 ",i,"th"
            elif i in lcid1:
                lc1=lc1 + allcts
            elif i in lcid2:
                lc2=lc2 + allcts
        else:
            if i in lcid0:
                lc0=np.append(lc0,allcts)
                #print lc0[0],lc0[1],lc0[2],"lc0 1st sect"
            elif i in lcid1:
                lc1=np.append(lc1,allcts)
            elif i in lcid2:
                lc2=np.append(lc2,allcts)
        hd.close()
        del hd
    hd0=pf.open(lcpath+infilekey+"_g%s_%s.lc"%(0,0))
    hd1=pf.open(lcpath+infilekey+"_g%s_%s.lc"%(0,1))
    hd2=pf.open(lcpath+infilekey+"_g%s_%s.lc"%(0,5))
    if (len(hd0[1].data.field(0))!=len(lc0)):
        print('together lc error!')
        sys.exit(0)
    col01 = pf.Column(name = 'Time', format='D',array = hd0[1].data.field(0))
    col02 = pf.Column(name = 'Counts', format='E',array = lc0)
    cols = pf.ColDefs([col01,col02])
    tbhdu = pf.TableHDU.from_columns(cols)
    tbhdu.writeto(lcpath+infilekey+"_lc0.fits",overwrite=True)###

    col11 = pf.Column(name = 'Time', format='D',array = hd1[1].data.field(0))
    col12 = pf.Column(name = 'Counts', format='E',array = lc1)
    cols = pf.ColDefs([col11,col12])
    tbhdu = pf.TableHDU.from_columns(cols)
    tbhdu.writeto(lcpath+infilekey+"_lc1.fits",overwrite=True)###

    col21 = pf.Column(name = 'Time', format='D',array = hd2[1].data.field(0))
    col22 = pf.Column(name = 'Counts', format='E',array = lc2)
    cols = pf.ColDefs([col21,col22])
    tbhdu = pf.TableHDU.from_columns(cols)
    tbhdu.writeto(lcpath+infilekey+"_lc5.fits",overwrite=True)###
    hd0.close()
    hd1.close()
    hd2.close()

def bkg_sub(x,y,m,flag=None):
    n = y.shape[0]
    v = np.zeros(n)
    w = np.zeros(n)
    b = np.zeros(n)
    v = np.log(np.log(np.sqrt(y+1)+1)+1)

    for p in range(m):
        for i in range(p+1,n-p):
            a1 = v[i]
            a2 = (v[i-p]+v[i+p])/2
            if flag==None:
                w[i] = a2 if a1>a2 else a1
            else:
                w[i] = a2
        for i in range(p+1,n-p):
            v[i] = w[i]
    b = (np.exp(np.exp(v)-1)-1)**2-1
    bkg = b
    src = y-bkg
    return bkg,src

'''
def smooth(inp,window_len=15,window='hanning'):
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    x = np.copy(inp)
    x_right = x[-(window_len-1+50):] - (np.sum(x[-60:-45])/15 - np.sum(x[-15:])/15)
    x_left = x[:window_len-1+50] - (np.sum(x[45:60])/15 - np.sum(x[0:15])/15)
    s=np.r_[x_left,x,x_right]
    if window == 'flat':
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    smtcut = (window_len-1)/2 + 50
    y=np.convolve(w/w.sum(),s,mode='valid')[smtcut:-smtcut]
    return y
'''
def smooth(inp,window_len=15,window='hanning'):
    x = np.copy(inp)
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<5:
        return x
    x_right = x[-(window_len-1):] - (np.sum(x[-window_len:-window_len+5])/5 - np.sum(x[-5:])/5)
    x_left = x[:window_len-1] - (np.sum(x[window_len-5:window_len])/5 - np.sum(x[0:5])/5)
    print('x_ritht,x_left',x_right,x_right.shape,x_left,x_left.shape)
    s=np.r_[x_left,x,x_right]
    if window == 'flat':
        w=np.ones(window_len,'d')
        print('if,',w,w.shape)
    else:
        w=eval('np.'+window+'(window_len)')
        print('else,', w, w.shape)
    smtcut = (window_len-1)/2
    print('smtcut:',smtcut)###
    y=np.convolve(w/w.sum(),s,mode='valid')[int(smtcut):-int(smtcut)]
    return y


def Smooth_Sysm(inp,window_len=11,window='hanning'):
    x = np.copy(inp)
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<5:
        return x
    x_right = x[-3:].mean()*2 - x[-1:-(window_len-1):-1]
    x_left = x[:3].mean()*2 - x[(window_len-1):1:-1]
    s=np.r_[x_left,x,x_right]
    if window == 'flat':
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    smtcut = (window_len-1)/2 - 1
    y=np.convolve(w/w.sum(),s,mode='valid')[int(smtcut):-int(smtcut)]
    return y

def enlarge(inp,window_len=11):
    x = np.copy(inp)
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<5:
        return x
    x_right = x[-3:].mean()*2 - x[-1:-(window_len-1):-1]
    x_left = x[:3].mean()*2 - x[(window_len-1):1:-1]
    s=np.r_[x_left,x,x_right]
    print(len(s))
    a=np.linspace(-1,1,window_len)
    w = np.cos(a*np.pi)
    w=w+(1.- w.sum())/float(window_len)
    print(w.sum())
    smtcut = int((window_len-1)/2-1)
    print(smtcut)
    y=np.convolve(w/w.sum(),s,mode='valid')[smtcut:-smtcut]
    print(len(y))
    return y

def Segments_Smooth(inplc,index,window_len=11):
    nub = len(index)
    smoothlc = np.array([])
    idx = np.array([])
    for i in range(nub):
        try:
            #print inplc[index[i]].shape
            smlc = Smooth_Sysm(inplc[index[i][30:-30]],window_len=window_len)
            #print smlc.shape
            smoothlc = np.append(smoothlc,smlc)
            idx = np.append(idx,index[i][30:-30])
            #print smlc.shape,idx.shape
        except:
            pass
    return smoothlc,idx.astype(np.int)

def smooth_LC(inpTM,inpLC,window_len=11,window='hanning'):
    x = np.copy(inpTM)
    y = np.copy(inpLC)
    x_right = x[-1]*2 - x[-1:-(window_len-1+50):-1]
    x_left = x[0]*2 - x[(window_len-1+50):1:-1]
    s=np.r_[x_left,x,x_right]
    if window == 'flat':
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    smtcut = int((window_len-1)/2 + 49)
    y=np.convolve(w/w.sum(),s,mode='valid')[smtcut:-smtcut]
    return y

def poly_bkg(lcfile):
    hdu = pf.open(lcfile)
    tb = hdu[1].data
    x = tb.field(0)
    y = tb.field(1)
    xr = x[y>0]
    yr = y[y>0]
    Ygti = ((xr[1:]-xr[:-1])!=1)
    c=np.nonzero(Ygti)[0]
    pot=np.r_[0,c+1,Ygti.shape[0]]
    par = []
    for i in range(len(pot)-1):
        try:
            yr1 = yr[pot[i]:pot[i+1]][2:-2]
            xr1 = xr[pot[i]:pot[i+1]][2:-2]
        except IndexError:
            raw_input("Indexerror ,recompose you code!")
            yr1 = yr[pot[i]:pot[i+1]]
            xr1 = xr[pot[i]:pot[i+1]]
        plt.plot(xr1,yr1,'r-',alpha=0.5)
        print('frag ',i,'len :',yr1.shape[0])
        if yr1.shape[0]>=600:
            #ysmooth = smooth(yr1,window_len=11)
            #yr1[(yr1-ysmooth)>60] = ysmooth[(yr1-ysmooth)>60]
            ysmooth = smooth(yr1,window_len=11,window='blackman')
            fraglen = 30
            startnum = 300/fraglen
            nmb = yr1.shape[0]/fraglen+2*startnum
            steplen = 15
            lcdic = {}
            #plt.plot(xr1,yr1,'r-',alpha=0.5)
            for idx in range(int(nmb)):
                stat = (idx-startnum) * fraglen  if (idx-startnum)*fraglen > 0 else 0
                stop = (idx-startnum) * fraglen + 600 if (idx-startnum) * fraglen + 600 < yr1.shape[0] else yr1.shape[0]
                if stat>=yr1.shape[0]-300:
                    stat=yr1.shape[0]-300
                    break
                stat = int(stat)
                stop = int(stop)
                frag = ysmooth[int(stat):int(stop)]
                fragmin = ysmooth[int(stat):int(stop)]
                tmmin = xr1[int(stat):int(stop)]
                polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
                for ratio in np.arange(1.,4.):
                    while 1:
                        fragres = polmin(tmmin-xr1[0]) - fragmin
                        if (fragres < -np.sqrt(fragmin)/ratio).sum()>0:
                            fragmin = fragmin[fragres!=fragres.min()]
                            tmmin = tmmin[fragres!=fragres.min()]
                            polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
                        else:
                            break
                    #plt.plot(tmmin,fragmin,'g+',ms=8)
                    tstmin,tspmin=FindGap(tmmin,24)
                    #print len(tmmin),len(tstmin)
                    for ix,it in enumerate(tstmin):
                        if it!=xr1[stat]:
                            tstmin[ix] = tstmin[ix]+10
                        if tspmin[ix] != xr1[stop-1]:
                            tspmin[ix] = tspmin[ix]-10
                    tmyep = GenArraFromGap(tstmin,tspmin).astype(np.int)
                    #print len(tmyep)
                    fragmin = fragmin[tmyep]
                    tmmin = tmmin[tmyep]
                    #print len(tmmin)
                    #plt.plot(tmmin,fragmin,'b.',ms=8)
                    #plt.show()
                    if len(tmmin)==0:continue
                    polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
                    while 1:
                        fragres = fragmin - polmin(tmmin-xr1[0])
                        if (fragres < -np.sqrt(fragmin)/ratio*1.1).sum()>0:
                            fragmin = fragmin[fragres!=fragres.min()]
                            tmmin = tmmin[fragres!=fragres.min()]
                            polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
                        else:
                            break
                #polmin = np.poly1d(np.polyfit(tmmin[30:-30]-xr1[0], fragmin[30:-30], 2))
                for ix in range(len(xr1[stat:stop])):
                    if (str(xr1[stat:stop][ix])) in lcdic:
                        lcdic[str(xr1[stat:stop][ix])].append(polmin(xr1[stat:stop][ix]-xr1[0]))
                    else:
                        lcdic[str(xr1[stat:stop][ix])]= [polmin(xr1[stat:stop][ix]-xr1[0])]
                #plt.plot(xr1[stat:stop],polmin(xr1[stat:stop]-xr1[0]),'g-',lw=1)
                plt.plot(tmmin,fragmin,'b.',ms=8)
            for key,vals in lcdic.items():
                lcdic[key]=np.median(vals)
                #plt.hist(vals,10)
                #plt.show()
            tptm = sorted(lcdic.keys())
            tplc = np.array([lcdic[tmid] for tmid in tptm])
            tptms = np.array([float(tmid) for tmid in tptm])
            ystp = Smooth_Sysm(tplc,window_len=151,window='hamming')
            ysmooth = Smooth_Sysm(yr1,window_len=11)
            for step in range(0+steplen,len(xr1),steplen):
                frag = ysmooth[step:step+steplen]
                tp = frag.min()
                timeloc = xr1[step:step+steplen][frag==tp][0]
                if np.array(ystp)[tptms==timeloc][0]>tp:
                    headcut = timeloc - xr1[0]
                    break
            #print len(xr1),len(yr1),"---"
            plt.plot(tptms,tplc,'b-',lw=1.5)
            #plt.plot(tptms,ystp,'r-',lw=1.5)
            plt.plot(xr1,ysmooth,'k--',lw=1.5)
            #plt.show()
            for step in range(len(xr1)-1-steplen,0,-steplen):
                #print len(ysmooth),step - steplen,step,"***"
                frag = ysmooth[step - steplen:step]
                tp = np.argmin(frag)
                timeloc = xr1[step - steplen:step][tp]
                if np.array(ystp)[tptms==timeloc][0]>frag[tp]:
                    endcut = timeloc - xr1[-1]
                    break
            #print headcut,endcut,">>>>>>>>>>>>>>>>>>"
            idxs = np.arange(len(yr1)).astype(np.int)
            bkg1 = ystp[int(headcut+60):int(endcut-60)]
            xr1 = tptms[int(headcut+60):int(endcut-60)]
            yr1 = yr1[int(headcut+60):int(endcut-60)]
            #plt.plot(xr1,bkg1,'.')
            print(len(xr1))
            k1 = (bkg1[2:]-bkg1[:-2])/2.
            sk1 = smooth(k1,301 if len(k1)>300 else len(k1)-1)
            k2 = sk1[2:]-sk1[:-2]
            sk2 = np.abs(smooth(k2,301 if len(k1)>300 else len(k2)-1))
            xr1 = xr1[2:-2]
            yr1 = yr1[2:-2]
            bkg1 = bkg1[2:-2]
            #print len(sk2),len(xr1)
            while 1:
                ink = np.argmax(sk2)
                if (sk2[ink]>0.001):
                    inkstart = ink-300 if ink>300 else 0
                    inkstop = ink+300 if ink<len(sk2) else 0
                    #print len(sk2),len(xr1)
                    yep = ~np.in1d(range(len(sk2)),range(inkstart,inkstop))
                    sk2 = sk2[yep]
                    xr1 = xr1[yep]
                    yr1 = yr1[yep]
                    bkg1 = bkg1[yep]
                else:
                    break
                #print len(sk2),len(xr1)
            plt.plot(xr1,bkg1,'y-')
            if i ==0:
                bkg = bkg1
                xr_sel = xr1
                yr_sel = yr1
            if i >0:
                try:
                    bkg = np.r_[bkg,bkg1]
                    xr_sel = np.r_[xr_sel,xr1]
                    yr_sel = np.r_[yr_sel,yr1]
                except UnboundLocalError:
                    bkg = bkg1
                    xr_sel = xr1
                    yr_sel = yr1
        ###plt.show()
        plt.close('all')    
    hdu.close()
    del hdu
    return xr_sel,yr_sel,bkg

def poly_bkg_array(x_org,y_org):
    x = x_org.copy()
    y = y_org.copy()
    xr = x[y>0]
    yr = y[y>0]
    Ygti = ((xr[1:]-xr[:-1])>1)
    c=np.nonzero(Ygti)[0]
    pot=np.r_[0,c+1,Ygti.shape[0]]
    par = []
    for i in range(len(pot)-1):
        try:
            yr1 = yr[pot[i]:pot[i+1]][2:-2]
            xr1 = xr[pot[i]:pot[i+1]][2:-2]
        except IndexError:
            raw_input("Indexerror ,recompose you code!")
            yr1 = yr[pot[i]:pot[i+1]]
            xr1 = xr[pot[i]:pot[i+1]]
        plt.plot(xr1,yr1,'r-',alpha=0.5)
        print('frag ',i,'len :',yr1.shape[0])
        if yr1.shape[0]>=600:
            #ysmooth = smooth(yr1,window_len=11)
            #yr1[(yr1-ysmooth)>60] = ysmooth[(yr1-ysmooth)>60]
            ysmooth = smooth(yr1,window_len=11,window='blackman')
            fraglen = 30
            startnum = 300/fraglen
            nmb = yr1.shape[0]/fraglen+2*startnum
            steplen = 15
            lcdic = {}
            #plt.plot(xr1,yr1,'r-',alpha=0.5)
            for idx in range(int(nmb)):
                stat = (idx-startnum) * fraglen  if (idx-startnum)*fraglen > 0 else 0
                stop = (idx-startnum) * fraglen + 600 if (idx-startnum) * fraglen + 600 < yr1.shape[0] else yr1.shape[0]
                if stat>=yr1.shape[0]-300:
                    stat=yr1.shape[0]-300
                    break
                frag = ysmooth[stat:stop]
                fragmin = ysmooth[stat:stop]
                tmmin = xr1[stat:stop]
                polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
                for ratio in np.arange(1.,4.):
                    while 1:
                        fragres = polmin(tmmin-xr1[0]) - fragmin
                        if (fragres < -np.sqrt(fragmin)/ratio).sum()>0:
                            fragmin = fragmin[fragres!=fragres.min()]
                            tmmin = tmmin[fragres!=fragres.min()]
                            polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
                        else:
                            break
                    #plt.plot(tmmin,fragmin,'g+',ms=8)
                    tstmin,tspmin=FindGap(tmmin,24)
                    #print len(tmmin),len(tstmin)
                    for ix,it in enumerate(tstmin):
                        if it!=xr1[stat]:
                            tstmin[ix] = tstmin[ix]+10
                        if tspmin[ix] != xr1[stop-1]:
                            tspmin[ix] = tspmin[ix]-10
                    tmyep = GenArraFromGap(tstmin,tspmin).astype(np.int)
                    #print len(tmyep)
                    fragmin = fragmin[tmyep]
                    tmmin = tmmin[tmyep]
                    #print len(tmmin)
                    #plt.plot(tmmin,fragmin,'b.',ms=8)
                    #plt.show()
                    if len(tmmin)==0:continue
                    polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
                    while 1:
                        fragres = fragmin - polmin(tmmin-xr1[0])
                        if (fragres < -np.sqrt(fragmin)/ratio*1.1).sum()>0:
                            fragmin = fragmin[fragres!=fragres.min()]
                            tmmin = tmmin[fragres!=fragres.min()]
                            polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
                        else:
                            break
                #polmin = np.poly1d(np.polyfit(tmmin[30:-30]-xr1[0], fragmin[30:-30], 2))
                for ix in range(len(xr1[stat:stop])):
                    if (str(xr1[stat:stop][ix])) in lcdic:
                        lcdic[str(xr1[stat:stop][ix])].append(polmin(xr1[stat:stop][ix]-xr1[0]))
                    else:
                        lcdic[str(xr1[stat:stop][ix])]= [polmin(xr1[stat:stop][ix]-xr1[0])]
                #plt.plot(xr1[stat:stop],polmin(xr1[stat:stop]-xr1[0]),'g-',lw=1)
                plt.plot(tmmin,fragmin,'b.',ms=8)
            for key,vals in lcdic.items():
                lcdic[key]=np.median(vals)
                #plt.hist(vals,10)
                #plt.show()
            tptm = sorted(lcdic.keys())
            tplc = np.array([lcdic[tmid] for tmid in tptm])
            tptms = np.array([float(tmid) for tmid in tptm])
            ystp = Smooth_Sysm(tplc,window_len=151,window='hamming')
            ysmooth = Smooth_Sysm(yr1,window_len=11)
            for step in range(0+steplen,len(xr1),steplen):
                frag = ysmooth[step:step+steplen]
                tp = frag.min()
                timeloc = xr1[step:step+steplen][frag==tp][0]
                if np.array(ystp)[tptms==timeloc][0]>tp:
                    headcut = timeloc - xr1[0]
                    break
            #print len(xr1),len(yr1),"---"
            plt.plot(tptms,tplc,'b-',lw=1.5)
            #plt.plot(tptms,ystp,'r-',lw=1.5)
            plt.plot(xr1,ysmooth,'k--',lw=1.5)
            #plt.show()
            for step in range(len(xr1)-1-steplen,0,-steplen):
                #print len(ysmooth),step - steplen,step,"***"
                frag = ysmooth[step - steplen:step]
                tp = np.argmin(frag)
                timeloc = xr1[step - steplen:step][tp]
                if np.array(ystp)[tptms==timeloc][0]>frag[tp]:
                    endcut = timeloc - xr1[-1]
                    break
            #print headcut,endcut,">>>>>>>>>>>>>>>>>>"
            idxs = np.arange(len(yr1)).astype(np.int)
            bkg1 = ystp[int(headcut+60):int(endcut-60)]
            xr1 = tptms[int(headcut+60):int(endcut-60)]
            yr1 = yr1[int(headcut+60):int(endcut-60)]
            #plt.plot(xr1,bkg1,'.')
            print(len(xr1))
            k1 = (bkg1[2:]-bkg1[:-2])/2.
            sk1 = smooth(k1,301 if len(k1)>300 else len(k1)-1)
            k2 = sk1[2:]-sk1[:-2]
            sk2 = np.abs(smooth(k2,301 if len(k1)>300 else len(k2)-1))
            xr1 = xr1[2:-2]
            yr1 = yr1[2:-2]
            bkg1 = bkg1[2:-2]
            #print len(sk2),len(xr1)
            while 1:
                ink = np.argmax(sk2)
                if (sk2[ink]>0.001):
                    inkstart = ink-300 if ink>300 else 0
                    inkstop = ink+300 if ink<len(sk2) else 0
                    #print len(sk2),len(xr1)
                    yep = ~np.in1d(range(len(sk2)),range(inkstart,inkstop))
                    sk2 = sk2[yep]
                    xr1 = xr1[yep]
                    yr1 = yr1[yep]
                    bkg1 = bkg1[yep]
                else:
                    break
                #print len(sk2),len(xr1)
            plt.plot(xr1,bkg1,'y-')
            if i ==0:
                bkg = bkg1
                xr_sel = xr1
                yr_sel = yr1
            if i >0:
                try:
                    bkg = np.r_[bkg,bkg1]
                    xr_sel = np.r_[xr_sel,xr1]
                    yr_sel = np.r_[yr_sel,yr1]
                except UnboundLocalError:
                    bkg = bkg1
                    xr_sel = xr1
                    yr_sel = yr1
        ###plt.show()
        plt.close('all')    
    hdu.close()
    del hdu
    return xr_sel,yr_sel,bkg


def snip_bkg(lcfile,m):
    hdu = pf.open(lcfile)
    tb = hdu[1].data
    x = tb.field(0)
    y = tb.field(1)
    y2 = y.copy()
    for i, yi in enumerate(y2):
        if i ==0:
            continue
        if i == y.shape[0]-2:
            break
        if yi == 0:
            y[i-1]=0
            y[i+1]=0
   
    xr = x[y>0]
    yr = y[y>0]
    Ygti = ((xr[1:]-xr[:-1])!=1)
    c=np.nonzero(Ygti)[0]
    pot=np.r_[0,c+1,Ygti.shape[0]]
    print(pot,"pot::",i)
    for i in range(len(pot)-1): 
        try:
            yr1 = yr[pot[i]:pot[i+1]][2:-2]
            xr1 = xr[pot[i]:pot[i+1]][2:-2]
        except IndexError:
            raw_input("Indexerror ,recompose you code!")
            yr1 = yr[pot[i]:pot[i+1]]
            xr1 = xr[pot[i]:pot[i+1]]
        #plt.plot(xr1,yr1,'r-',alpha=0.5)
        if yr1.shape[0]>=301:
            print((yr1.shape),"len first yr1::",i)
            ysmooth = smooth(yr1,window_len=11)
            yr1[(yr1-ysmooth)>60] = ysmooth[(yr1-ysmooth)>60]
            ysmooth = smooth(yr1,window_len=11)
            #plt.plot(xr1,ysmooth,'k-',alpha=0.5)
            print((xr1.shape),(ysmooth.shape),"input to sub")
            bkg1,src1 = bkg_sub(xr1,ysmooth,m)
            #plt.plot(xr1,bkg1,'y-',alpha=1)
            bkg1,src1 = bkg_sub(xr1,bkg1,8,flag=1)
            loc1,loc2=0,0
            head = (ysmooth-bkg1)[:150]>15
            head = np.nonzero(head)[0]
            headsum = ((head[1:]-head[:-1])<3).sum()
            rate1 = (bkg1[1:61] - bkg1[:60])
            tail = (ysmooth-bkg1)[-150:]>15
            tail = np.nonzero(tail)[0]
            tailsum = ((tail[1:]-tail[:-1])<3).sum()
            rate2 = (bkg1[-60:]-bkg1[-61:-1])
            print("Test",headsum,tailsum)
            headsum = ((ysmooth-bkg1))[0:60].sum()
            tailsum = ((ysmooth-bkg1))[-60:].sum()
            print("Test",headsum,tailsum)
            if headsum > 600:
                for clen in range(300):
                    print("ail",clen,((ysmooth-bkg1)[(30+clen):(90+clen)]).sum())
                    if ((ysmooth-bkg1)[30+clen:90+clen]).sum()<500:
                        print('Head',clen,((ysmooth-bkg1)[30+clen:90+clen]).sum())
                        break
                locnum=np.nonzero(ysmooth[30+clen:60+clen]==ysmooth[30+clen:60+clen].min())[0]
                loc1 = 30+clen+locnum[0]
                print('locnum',locnum,ysmooth[30+clen:60+clen].min(),ysmooth[loc1])
            if tailsum > 600:
                for clen in range(300):
                    print("Tail",clen,((ysmooth-bkg1)[-(90+clen):-(30+clen)]).sum())
                    if ((ysmooth-bkg1)[-(90+clen):-(30+clen)]).sum()<500:
                        print("Tail",clen,((ysmooth-bkg1)[-(90+clen):-(30+clen)]).sum())
                        break
                locnum=np.nonzero(ysmooth[-(60+clen):-(30+clen)]==ysmooth[-(60+clen):-(30+clen)].min())[0]
                loc2 = 60+clen-locnum[0]
            if loc1+loc2>0:
                xr1=xr1[loc1:]
                yr1=yr1[loc1:]
                if loc2>0:
                    xr1 = xr1[:-loc2]
                    yr1 = yr1[:-loc2]
                ysmooth = smooth(yr1,window_len=11)
                yr1[(yr1-ysmooth)>60] = ysmooth[(yr1-ysmooth)>60]
                ysmooth = smooth(yr1,window_len=11)
                bkg1,src1 = bkg_sub(xr1,ysmooth,m)
                bkg1,src1 = bkg_sub(xr1,bkg1,8,flag=1)
            print("Tesr End!")
            cutpri = 90 if abs((bkg1[149]-bkg1[0])/150.0)< 0.5 else 149
            cutend = 90 if abs((bkg1[-150]-bkg1[-1])/150.0) < 0.5 else 149
            print(cutpri,(bkg1[149]-bkg1[0])/150.0,cutend,(bkg1[-90]-bkg1[-1])/90.0)
            test=(yr1-bkg1)
            #plt.plot(xr1,yr1-bkg1,'bo',alpha=0.3,ms=7)
            yth = smooth(yr1-bkg1,window_len=31)
            #plt.plot(xr1,yth,'r-',alpha=1)
            #plt.plot(xr1,bkg1,'g-',alpha=0.5)
            #plt.axhline(test.std()+test.mean(),ls='--',color='r',lw=1)
            #plt.axhline(test.mean()-test.std(),ls='--',color='r',lw=1)
            #plt.grid()
            #plt.show()
            plt.close('all')
            if i ==0:
                bkg = bkg1[cutpri:-cutend]
                xr_sel = xr1[cutpri:-cutend]
                yr_sel = yr1[cutpri:-cutend]
            if i >0:
                try:
                    bkg = np.r_[bkg,bkg1[cutpri:-cutend]]
                    xr_sel = np.r_[xr_sel,xr1[cutpri:-cutend]]
                    yr_sel = np.r_[yr_sel,yr1[cutpri:-cutend]]
                except UnboundLocalError:
                    bkg = bkg1[cutpri:-cutend]
                    xr_sel = xr1[cutpri:-cutend]
                    yr_sel = yr1[cutpri:-cutend]
        #print i,yr1.shape,bkg.shape
        if yr1.shape[0]<300:
            print("the gti is to short: ")
    #print "the shape of bkg:",bkg,bkg.shape
    try:
        xr_sel+0
    except NameError:
        print("There's No good times!")
        return 0
    return xr_sel,yr_sel,bkg

def snip(key):
    lcpath = "/hxmt/work/HXMT_scan_data/HE/nydata2/lc/" + key
    bkgpath = "/hxmt/work/HXMT_scan_data/HE/nydata2/bkg/"
    pst = "snip"
    print(lcpath)

    xr0_small,yr0_small,bkg0_small = snip_bkg(lcpath+"_lc0.fits",49)
    print("me 0 done")
    xr1_small,yr1_small,bkg1_small = snip_bkg(lcpath+"_lc1.fits",47)
    print("me 1 done")
    xr2_small,yr2_small,bkg2_small = snip_bkg(lcpath+"_lc5.fits",49)
    unitime = np.array(Uniform([xr0_small.tolist(),xr1_small.tolist(),xr2_small.tolist()]))
    hd = pf.open(lcpath[:-4]+"Ndead_lc0.fits")
    fitstm0 = hd[1].data.field(0);fitrate0 = hd[1].data.field(1)
    hd.close();del hd
    hd = pf.open(lcpath[:-4]+"Ndead_lc1.fits")
    fitstm1 = hd[1].data.field(0);fitrate1 = hd[1].data.field(1)
    hd.close();del hd
    hd = pf.open(lcpath[:-4]+"Ndead_lc5.fits")
    fitstm2 = hd[1].data.field(0);fitrate2 = hd[1].data.field(1)
    hd.close();del hd
    rat0 = fitrate0[np.in1d(fitstm0,unitime)]
    rat1 = fitrate1[np.in1d(fitstm1,unitime)]
    rat2 = fitrate2[np.in1d(fitstm2,unitime)]
    yep0 = np.in1d(xr0_small,unitime)
    yep1 = np.in1d(xr1_small,unitime)
    yep2 = np.in1d(xr2_small,unitime)

    yr0 = yr0_small[yep0]
    bkg0 = bkg0_small[yep0]
    err0 = np.sqrt(rat0)*(yr0/rat0)
    col1 = pf.Column(name = 'Time', format='D',array = unitime)
    col2 = pf.Column(name = 'Counts', format='D',array = yr0-bkg0)
    col3 = pf.Column(name = 'Stat_err', format='D',array = err0)
    col4 = pf.Column(name = 'bkg', format='D',array = bkg0)
    cols = pf.ColDefs([col1,col2,col3,col4])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    tbhdu.writeto(bkgpath+key + pst + "0.fits",clobber=True)

    yr1 = yr1_small[yep1]
    bkg1 = bkg1_small[yep1]
    err1 = np.sqrt(rat1)*(yr1/rat1)
    col1 = pf.Column(name = 'Time', format='D',array = unitime)
    col2 = pf.Column(name = 'Counts', format='D',array = yr1-bkg1)
    col3 = pf.Column(name = 'Stat_err', format='D',array = err1)
    col4 = pf.Column(name = 'bkg', format='D',array = bkg1)
    cols = pf.ColDefs([col1,col2,col3,col4])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    tbhdu.writeto(bkgpath+key + pst + "1.fits",clobber=True)

    yr2 = yr2_small[yep2]
    bkg2 = bkg2_small[yep2]
    err2 = np.sqrt(rat2)*(yr2/rat2)
    col1 = pf.Column(name = 'Time', format='D',array = unitime)
    col2 = pf.Column(name = 'Counts', format='D',array = yr2-bkg2)
    col3 = pf.Column(name = 'Stat_err', format='D',array = err2)
    col4 = pf.Column(name = 'bkg', format='D',array = bkg2)
    cols = pf.ColDefs([col1,col2,col3,col4])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    tbhdu.writeto(bkgpath+key + pst + "5.fits",clobber=True)

    fig,axe = plt.subplots(3,sharex=True, sharey=False,figsize=(22,16))
    fig.subplots_adjust(hspace = 0.05,top=0.9,bottom=0.1,left=0.1,right=0.95)
    l0=axe[0].errorbar(unitime,yr0,yerr=err0,fmt='.',c='k',ecolor='b')
    l1=axe[1].errorbar(unitime,yr1,yerr=err1,fmt='.',c='k',ecolor='b')
    l2=axe[2].errorbar(unitime,yr2,yerr=err2,fmt='.',c='k',ecolor='b')
    l3=axe[0].plot(unitime,bkg0,c='r',lw=1.5, zorder=50)
    l4=axe[1].plot(unitime,bkg1,c='r',lw=1.5, zorder=50)
    l5=axe[2].plot(unitime,bkg2,c='r',lw=1.5, zorder=50)
    l6=axe[0].plot(unitime,yr0-bkg0,c='k',lw=0.5)
    l7=axe[1].plot(unitime,yr1-bkg1,c='k',lw=0.5)
    l8=axe[2].plot(unitime,yr2-bkg2,c='k',lw=0.5)
    axe[1].set_ylabel('counts/sec',fontsize=15)

    for i in range(3):
        axe[i].grid()
        axe[i].legend(loc='upper right')
        axe[i].legend(labels=['data','background','residual'], loc='upper left')

    axe[0].set_title(key[:13]+"_he")
    plt.xlabel('Time /s',fontsize = 15)
    plt.savefig(bkgpath+key+"snip_bkg.png",dpi=400)
    #plt.show()
    plt.close('all')

def bkg_muti_poly(key,inpath,bkgpath,pst = '_netlc'):
    lcpath = inpath + key
    xr0_small,yr0_small,bkg0_small = poly_bkg(lcpath+"_lc0.fits");print("me 0 done")###
    xr1_small,yr1_small,bkg1_small = poly_bkg(lcpath+"_lc1.fits");print("me 1 done")###
    xr2_small,yr2_small,bkg2_small = poly_bkg(lcpath+"_lc5.fits");print("me 2 done")###
    unitime = np.array(Uniform([xr0_small.tolist(),xr1_small.tolist(),xr2_small.tolist()]))
    hd = pf.open(lcpath[:-4]+"Ndead_lc0.fits")
    fitstm0 = hd[1].data.field(0);fitrate0 = hd[1].data.field(1)
    hd.close();del hd
    hd = pf.open(lcpath[:-4]+"Ndead_lc1.fits")
    fitstm1 = hd[1].data.field(0);fitrate1 = hd[1].data.field(1)
    hd.close();del hd
    hd = pf.open(lcpath[:-4]+"Ndead_lc5.fits")
    fitstm2 = hd[1].data.field(0);fitrate2 = hd[1].data.field(1)
    hd.close();del hd
    rat0 = fitrate0[np.in1d(fitstm0,unitime)]
    rat1 = fitrate1[np.in1d(fitstm1,unitime)]
    rat2 = fitrate2[np.in1d(fitstm2,unitime)]
    yep0 = np.in1d(xr0_small,unitime)
    yep1 = np.in1d(xr1_small,unitime)
    yep2 = np.in1d(xr2_small,unitime)

    yr0 = yr0_small[yep0]
    bkg0 = bkg0_small[yep0]
    err0 = np.sqrt(rat0)*(yr0/rat0)
    col1 = pf.Column(name = 'Time', format='D',array = unitime)
    col2 = pf.Column(name = 'Counts', format='D',array = yr0-bkg0)
    col3 = pf.Column(name = 'Stat_err', format='D',array = err0)
    col4 = pf.Column(name = 'bkg', format='D',array = bkg0)
    cols = pf.ColDefs([col1,col2,col3,col4])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    tbhdu.writeto(bkgpath+key + pst + "_b0.fits",clobber=True)

    yr1 = yr1_small[yep1]
    bkg1 = bkg1_small[yep1]
    err1 = np.sqrt(rat1)*(yr1/rat1)
    col1 = pf.Column(name = 'Time', format='D',array = unitime)
    col2 = pf.Column(name = 'Counts', format='D',array = yr1-bkg1)
    col3 = pf.Column(name = 'Stat_err', format='D',array = err1)
    col4 = pf.Column(name = 'bkg', format='D',array = bkg1)
    cols = pf.ColDefs([col1,col2,col3,col4])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    tbhdu.writeto(bkgpath+key + pst + "_b1.fits",clobber=True)###

    yr2 = yr2_small[yep2]
    bkg2 = bkg2_small[yep2]
    err2 = np.sqrt(rat2)*(yr2/rat2)
    col1 = pf.Column(name = 'Time', format='D',array = unitime)
    col2 = pf.Column(name = 'Counts', format='D',array = yr2-bkg2)
    col3 = pf.Column(name = 'Stat_err', format='D',array = err2)
    col4 = pf.Column(name = 'bkg', format='D',array = bkg2)
    cols = pf.ColDefs([col1,col2,col3,col4])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    tbhdu.writeto(bkgpath+key + pst + "_b2.fits",clobber=True)###

    fig,axe = plt.subplots(3,sharex=True, sharey=False,figsize=(22,16))
    fig.subplots_adjust(hspace = 0.05,top=0.9,bottom=0.1,left=0.1,right=0.95)
    axe[0].errorbar(fitstm0,fitrate0,yerr=np.sqrt(fitrate0),fmt='.',c='k',ecolor='r',ms=5)
    axe[1].errorbar(fitstm1,fitrate1,yerr=np.sqrt(fitrate1),fmt='.',c='k',ecolor='r',ms=5)
    axe[2].errorbar(fitstm2,fitrate2,yerr=np.sqrt(fitrate2),fmt='.',c='k',ecolor='r',ms=5)

    l0=axe[0].errorbar(unitime,yr0,yerr=err0,fmt='.',c='k',ecolor='b',ms=5)
    l1=axe[1].errorbar(unitime,yr1,yerr=err1,fmt='.',c='k',ecolor='b',ms=5)
    l2=axe[2].errorbar(unitime,yr2,yerr=err2,fmt='.',c='k',ecolor='b',ms=5)
    l6=axe[0].plot(unitime,yr0-bkg0,c='k',lw=0.5)
    l7=axe[1].plot(unitime,yr1-bkg1,c='k',lw=0.5)
    l8=axe[2].plot(unitime,yr2-bkg2,c='k',lw=0.5)
    l3=axe[0].plot(unitime,bkg0,c='r',lw=1.5, zorder=50)
    l4=axe[1].plot(unitime,bkg1,c='r',lw=1.5, zorder=50)
    l5=axe[2].plot(unitime,bkg2,c='r',lw=1.5, zorder=50)
    axe[1].set_ylabel('counts/sec',fontsize=15)

    for i in range(3):
        axe[i].grid()
        #axe[i].legend(loc='upper right')
        axe[i].legend(labels=['net counts','background','N-Dead','data'], loc='upper left',fontsize=10)

    axe[0].set_title(key[:13]+"_he")
    plt.xlabel('Time /s',fontsize = 15)
    plt.savefig(bkgpath+key+"small_bkg.png",dpi=400)
    #plt.show()
    plt.close('all')

def correct(filename,outfile = './crtatt.fits'):
    hd = pf.open(filename)
    tm1 = hd[1].data.field(0)
    jd =np.array([0.0,0.25,0.50,0.75,1.0])
    mx0 = (abs(tm1 - np.floor(tm1)) <= 0.125)
    mx1 = (abs(tm1 - np.floor(tm1) - jd[1]) <= 0.125)
    mx2 = (abs(tm1 - np.floor(tm1) - jd[2]) <= 0.125) 
    mx3 = (abs(tm1 - np.floor(tm1) - jd[3]) <= 0.125)
    mx4 = (abs(tm1 - np.floor(tm1) - jd[4]) <= 0.125)
    tm1[mx0] = np.floor(tm1[mx0])
    tm1[mx1] = np.floor(tm1[mx1]) + 0.25
    tm1[mx2] = np.floor(tm1[mx2]) + 0.50
    tm1[mx3] = np.floor(tm1[mx3]) + 0.75
    tm1[mx4] = np.floor(tm1[mx4]) + 1.0
    b = ((tm1[1:]-tm1[:-1])>1)
    c = np.nonzero(b)[0]
    potstart = np.r_[0,c+1]
    potstop = np.r_[c,tm1.shape[0]-1]
    start= tm1[potstart]
    stop = tm1[potstop]
    temp = np.array([])
    for i in range(len(stop)):
        temp = np.append(temp,np.arange(start[i],stop[i]+0.2,0.25))
    tmt = temp.copy()
    for i in range(1,len(hd),1):
        temp=[]
        tb = hd[i].data
        temp.append(list(tmt))
        tmold = np.array((tb.field(0)).tolist())
        for j in range(1,len(tb.formats),1):
            #print i,j
            data = tb.field(j)
            newdata=np.interp(tmt,tmold,data)
            temp.append(list(newdata))
        tplt = list(zip(*temp))###np.array(zip(*temp)).tolist()
        fitdata = np.rec.array(tplt,names=','.join(tb.names))
        hd[i].data = fitdata
    hd.writeto(outfile,clobber=True)
    hd.close()
    del hd

def plotnetlc(filenames,outfile):
    key = filenames[0].split('/')[-1][:13]
    lchd0 = pf.open(filenames[0])
    lchd1 = pf.open(filenames[1])
    lchd2 = pf.open(filenames[2])
    
    lctb0 = lchd0[1].data
    unitime = lctb0['Time']
    net0 = lctb0['Counts']
    net1 = lchd1[1].data['Counts']
    net2 = lchd2[1].data['Counts']
    bkg0 = lctb0['bkg']
    bkg1 = lchd1[1].data['bkg']
    bkg2 = lchd2[1].data['bkg']
    err0 = lctb0['Stat_err']
    err1 = lchd1[1].data['Stat_err']
    err2 = lchd2[1].data['Stat_err']
    yr0 = net0 + bkg0
    yr1 = net1 + bkg1
    yr2 = net2 + bkg2
    fig,axe = plt.subplots(3,sharex=True, sharey=False,figsize=(22,16))
    fig.subplots_adjust(hspace = 0.05,top=0.9,bottom=0.1,left=0.1,right=0.95)
    l0=axe[0].errorbar(unitime,yr0,yerr=err0,fmt='.',c='k',ecolor='b',ms=5)
    l1=axe[1].errorbar(unitime,yr1,yerr=err1,fmt='.',c='k',ecolor='b',ms=5)
    l2=axe[2].errorbar(unitime,yr2,yerr=err2,fmt='.',c='k',ecolor='b',ms=5)
    l6=axe[0].plot(unitime,net0,c='k',lw=0.5)
    l7=axe[1].plot(unitime,net1,c='k',lw=0.5)
    l8=axe[2].plot(unitime,net2,c='k',lw=0.5)
    l3=axe[0].plot(unitime,bkg0,c='r',lw=1.5, zorder=50)
    l4=axe[1].plot(unitime,bkg1,c='r',lw=1.5, zorder=50)
    l5=axe[2].plot(unitime,bkg2,c='r',lw=1.5, zorder=50)
    axe[1].set_ylabel('counts  /sec',fontsize=15)
    
    for i in range(3):
        axe[i].grid()
        #axe[i].legend(loc='upper right')
        axe[i].legend(labels=['net counts','background','data'], loc='upper left',fontsize=10)
    
    axe[0].set_title(key[:13]+"_he")
    plt.xlabel('Time /s',fontsize = 15)
    plt.savefig(outfile,dpi=400)
    #plt.show()
    plt.close('all')
    lchd0.close()
    lchd1.close()
    lchd2.close()

