#!/hxmt/home/lick/soft/anaconda2/bin/python
'''
Copyright 2016, HXMT Ground System
All right reserved
File Name: hxmtpsf.py for sample
.......
Author: lick
Version: 0.1
Date:  2016-12-11
Email: lick@ihep.ac.cn
......
'''
from astropy.io import fits as pf
import numpy as np
from math import *
from readxml import *
import time
import Quat as quat
import matplotlib.pyplot as plt
import math
#############read ra dec for scan stripe################
def src_map(cfg):
#if len(sys.argv)<2:
#    print "Need the config file!"
#    sys.exit(1)
#print sys.argv[0],sys.argv[1]
#cfg = sys.argv[1]
    readcfg = loadDom(cfg)
    infilestr = readcfg.getTagText("infilenamelist")
    inpathstr = readcfg.getTagText("inpath")
    outpathstr = readcfg.getTagText("outpath")
    outpathstr = outpathstr.strip()
    outfilestr0 = outpathstr.split()[0]#.encode()
    #print inpathstr,infilestr
    inpathstr = inpathstr.strip()
    evtfilestr = infilestr.split()[0]
    evtfilestr2 = infilestr.split()[-2]
    infilestr2 = infilestr.split()[-1]
    
    infile = (inpathstr+infilestr2)#.encode()
    evtfile = (inpathstr+evtfilestr)#.encode()
    evtfile2 = (inpathstr+evtfilestr2)#.encode()
    filekey = infile.split('/')[-2][0:13]
    outfilestr0 = outfilestr0 +"/"+ filekey
    print ("the scanning pointing file : ",infile,evtfile)
    print ("the LC file : ",evtfile)
    hdulist = pf.open(infile)
    evtshdu = pf.open(evtfile)
    evtshdu2 = pf.open(evtfile2)
    tb1 = hdulist[1].data
    tb3 = hdulist[3].data
    tb2 = evtshdu[1].data
    tb4 = evtshdu2[1].data
    qtime = tb3.field(0)
    etime1 = tb2.field(0)
    etime2 = tb4.field(0)
    tstart = etime1[0]
    tstop = etime2[-1]
    #etime = np.r_[etime1,etime2]
    etime = etime1
    #etime = etime[etime>173343100]
    print("tstart,tstop: ",tstart,tstop)
    #tb3sel = tb3[(tb3.field(0)<=tstop)&(tb3.field(0)>=tstart)]
    tb3sel = tb3[(tb3.field(0)<=tstop) & (tb3.field(0)>=tstart)]
    print(tb3sel.shape)
    time.sleep(5)
    q1_list = tb3.field(1)[np.in1d(qtime,etime)]
    q2_list = tb3.field(2)[np.in1d(qtime,etime)]
    q3_list = tb3.field(3)[np.in1d(qtime,etime)]
    print(qtime[0:10],etime[0:10])
    stripe_ra_list = []
    stripe_dec_list = []
    print("q list size: ",q1_list.shape,sum(np.in1d(qtime,etime)))
    fatt=open(outfilestr0+"att.dat",'w')
    fq=open(outfilestr0+"q.dat",'w')
    for i in range(0,q1_list.shape[0]):
        qu = [q1_list[i],q2_list[i],q3_list[i]]
        qu = qu+[sqrt(fabs(1-qu[0]*qu[0]-qu[1]*qu[1]-qu[2]*qu[2]))]
        quat1 = quat.Quat(quat.normalize(qu))
        stripe_ra_list.append(quat1.ra)
        stripe_dec_list.append(quat1.dec)
        #print "ra dec: "
        #print quat1.ra,quat1.dec,quat1.roll
        print(quat1.ra,quat1.dec,quat1.roll,file=fatt)
        print(q1_list[i],q2_list[i],q3_list[i],file=fq)
    time.sleep(5)
    fatt.close()
    fq.close()
    #stripe_ra_list = t1.field(1)[(qtime<=tstop) & (qtime>=tstart)]
    #stripe_dec_list = t1.field(2)[(qtime<=tstop) & (qtime>=tstart)]
    stripe_ra_list = np.array(stripe_ra_list)
    stripe_dec_list = np.array(stripe_dec_list)
    ra_centre = math.atan2(np.sin(stripe_ra_list * math.pi/180.).mean(),np.cos(stripe_ra_list * math.pi/180.).mean())*180./math.pi
    ra_centre = ra_centre + 360. if ra_centre<0 else ra_centre
    dec_centre = stripe_dec_list.mean()
    #ra_centre = 270.
    #dec_centre = -30.0
    ralb = stripe_ra_list.min()
    rahb = stripe_ra_list.max()
    declb = stripe_dec_list.min()
    dechb = stripe_dec_list.max()
    
    src_list = pf.open("/hxmt/work/HXMT_scan_data/HE/lick_tools/gnrl_refr_cat_0031.fits")
    src_map = src_list[1].data
    src_name = src_map.field(2)
    src_ra = src_map.field(4)
    src_dec = src_map.field(5)
    src_flux = src_map.field(28)
    Yary=[]
    for i,itra in enumerate(src_ra):
        #print ((itra - stripe_ra_list)**2<34).sum(),(itra - stripe_ra_list)**2
        Yseq = (np.sin(src_dec[i]*np.pi/360.0 - stripe_dec_list*np.pi/360.0))**2 + np.cos(src_dec[i]*np.pi/180.0)*np.cos(stripe_dec_list*np.pi/180.0)*(np.sin(itra*np.pi/360.0 - stripe_ra_list*np.pi/360.0))**2
        Yseq = np.arccos(1-2*Yseq)*180./np.pi
        Yseq = (Yseq<5.7).sum()
        if Yseq>0:
            Yary.append(1)
        else:
            Yary.append(0)
    #Yary=((src_ra[:]<=rahb+6)&(src_ra[:]>=ralb-6)&(src_dec[:]>=declb-6)&(src_dec[:]<=dechb+6))
    #print sum(Yary)
    Yary=(np.array(Yary)>0)
    src_num = np.nonzero(Yary)[0]
    src_ra_sel = src_ra[Yary]
    src_dec_sel = src_dec[Yary]
    src_flux_sel = src_flux[Yary]
    src_name_sel = src_name[Yary]
    #src_flux_tot = np.sum(src_flux_sel,axis=1) ###if use the 32 colunms of flux 
    src_flux_tot = src_flux_sel
    #src_map_sel = np.vstack((src_flux_tot,src_ra_sel,src_dec_sel,src_name_sel))
    src_map_sel = np.vstack((src_flux_tot,src_ra_sel,src_dec_sel))
    print(src_map_sel.shape,len(src_ra_sel),len(src_dec_sel),len(src_flux_sel),len(src_name_sel))
    np.savetxt(outfilestr0+"integral.dat",src_map_sel.T,fmt="%s")
    return src_map_sel.T,ralb,rahb,declb,dechb,ra_centre,dec_centre,src_name_sel,src_num
    #for i, src_ra in enumerate(src_ra):

def src_new(cfg):
    readcfg = loadDom(cfg)
    infilestr = readcfg.getTagText("infilenamelist")
    inpathstr = readcfg.getTagText("inpath")
    outpathstr = readcfg.getTagText("outpath")
    instr = readcfg.getTagText("Instrument")
    outpathstr = outpathstr.strip()
    outfilestr0 = outpathstr.split()[0]#.encode()
    #print inpathstr,infilestr
    inpathstr = inpathstr.strip()
    evtfilestr = infilestr.split()[0]
    evtfilestr2 = infilestr.split()[-2]
    infilestr2 = infilestr.split()[-1]
     
    infile = (inpathstr+infilestr2)#.encode()
    evtfile = (inpathstr+evtfilestr)#.encode()
    evtfile2 = (inpathstr+evtfilestr2)#.encode()
    filekey = infile.split('/')[-2][0:13]
    outfilestr0 = outfilestr0 +"/"+ filekey
    print ("the scanning pointing file : ",infile,evtfile)
    print ("the LC file : ",evtfile)
    instr = instr.strip()
    instr = instr.split()[0]
    instr = instr#.encode()
    if instr == "HE":
        big_range=5.7 #means long_bound
        small_range = 1.1 #means short_bound
    elif instr == "ME":
        big_range=5.7 #means long_bound
        small_range = 1.1 #means short_bound
    elif instr == "LE":
        big_range=5.7 #means long_bound
        small_range = 1.1 #means short_bound
    else:
        print("The instrment in config file is wrong.")
        sys.exit(1)

    hdulist = pf.open(infile)
    evtshdu = pf.open(evtfile)
    evtshdu2 = pf.open(evtfile2)
    tb1 = hdulist[1].data
    tb3 = hdulist[3].data
    tb2 = evtshdu[1].data
    tb4 = evtshdu2[1].data
    qtime = tb3.field(0)
    etime1 = tb2.field(0)
    etime2 = tb4.field(0)
    tstart = etime1[0]
    tstop = etime2[-1]
    #etime = np.r_[etime1,etime2]
    etime = etime1
    #etime = etime[etime>173343100]
    print("tstart,tstop: ",tstart,tstop)
    #tb3sel = tb3[(tb3.field(0)<=tstop)&(tb3.field(0)>=tstart)]
    tb3sel = tb3[(tb3.field(0)<=tstop) & (tb3.field(0)>=tstart)]
    print(tb3sel.shape)
    time.sleep(5)
    q1_list = tb3.field(1)[np.in1d(qtime,etime)]
    q2_list = tb3.field(2)[np.in1d(qtime,etime)]
    q3_list = tb3.field(3)[np.in1d(qtime,etime)]
    qtime_uni = qtime[np.in1d(qtime,etime)]
    stripe_ra_list = []
    stripe_dec_list = []
    print("q list size: ",q1_list.shape,sum(np.in1d(qtime,etime)))
    fatt=open(outfilestr0+"att.dat",'w')
    fq=open(outfilestr0+"q.dat",'w')
    for i in range(0,q1_list.shape[0]):
        qu = [q1_list[i],q2_list[i],q3_list[i]]
        qu = qu+[sqrt(fabs(1-qu[0]*qu[0]-qu[1]*qu[1]-qu[2]*qu[2]))]
        quat1 = quat.Quat(quat.normalize(qu))
        stripe_ra_list.append(quat1.ra)
        stripe_dec_list.append(quat1.dec)
        #print "ra dec: "
        #print quat1.ra,quat1.dec,quat1.roll
        print(quat1.ra,quat1.dec,quat1.roll,fatt)
        print(q1_list[i],q2_list[i],q3_list[i],fq)
    time.sleep(5)
    fatt.close()
    fq.close()
    #stripe_ra_list = t1.field(1)[(qtime<=tstop) & (qtime>=tstart)]
    #stripe_dec_list = t1.field(2)[(qtime<=tstop) & (qtime>=tstart)]
    stripe_ra_list = np.array(stripe_ra_list)
    stripe_dec_list = np.array(stripe_dec_list)
    np.savetxt('sate_Info.txt',[qtime_uni.tolist(),stripe_ra_list.tolist(),stripe_dec_list.tolist(),q1_list.tolist(),q2_list.tolist(),q3_list.tolist()])
    ra_centre = math.atan2(np.sin(stripe_ra_list * math.pi/180.).mean(),np.cos(stripe_ra_list * math.pi/180.).mean())*180./math.pi
    ra_centre = ra_centre + 360. if ra_centre<0 else ra_centre
    dec_centre = stripe_dec_list.mean()
    ralb = stripe_ra_list.min()
    rahb = stripe_ra_list.max()
    declb = stripe_dec_list.min()
    dechb = stripe_dec_list.max()
    
    src_list = pf.open("/sharefs/hbkg/user/nangyi/lick_tools/para_files/Srcs_IGR_SWIFT_06.fits")
    src_map = src_list[1].data
    src_name = src_map.field('NAME')
    src_ra = src_map.field('RA')
    src_dec = src_map.field('DEC')
    src_flux = src_map.field('LE_NORM')
    Yary=[]
    for i,itra in enumerate(src_ra):
        #print ((itra - stripe_ra_list)**2<34).sum(),(itra - stripe_ra_list)**2
        Yseq = (np.sin(src_dec[i]*np.pi/360.0 - stripe_dec_list*np.pi/360.0))**2 + np.cos(src_dec[i]*np.pi/180.0)*np.cos(stripe_dec_list*np.pi/180.0)*(np.sin(itra*np.pi/360.0 - stripe_ra_list*np.pi/360.0))**2
        Yseq = np.arccos(1-2*Yseq)*180./np.pi
        Yseq = (Yseq<5.7).sum()
        if Yseq>0:
            Yary.append(1)
        else:
            Yary.append(0)
    print(sum(Yary))
    Yary=(np.array(Yary)>0)
    src_num = np.nonzero(Yary)[0]
    src_ra_sel = src_ra[Yary]
    src_dec_sel = src_dec[Yary]
    src_flux_sel = src_flux[Yary]
    src_name_sel = src_name[Yary]
    #src_flux_tot = np.sum(src_flux_sel,axis=1) ###if use the 32 colunms of flux 
    src_flux_tot = src_flux_sel
    #src_map_sel = np.vstack((src_flux_tot,src_ra_sel,src_dec_sel,src_name_sel))
    src_map_sel = np.vstack((src_flux_tot,src_ra_sel,src_dec_sel))
    print(src_map_sel.shape,len(src_ra_sel),len(src_dec_sel),len(src_flux_sel),len(src_name_sel))
    Yary=[]
    for i,itra in enumerate(src_ra):
        #print ((itra - stripe_ra_list)**2<34).sum(),(itra - stripe_ra_list)**2
        Yseq = (np.sin(src_dec[i]*np.pi/360.0 - stripe_dec_list*np.pi/360.0))**2 + np.cos(src_dec[i]*np.pi/180.0)*np.cos(stripe_dec_list*np.pi/180.0)*(np.sin(itra*np.pi/360.0 - stripe_ra_list*np.pi/360.0))**2
        Yseq = np.arccos(1-2*Yseq)*180./np.pi
        Yseq = (Yseq<5.81).sum()
        if Yseq>0:
            Yary.append(1)
        else:
            Yary.append(0)
    print(sum(Yary))
    Yary=(np.array(Yary)>0)
    src_ra_near = src_ra[Yary]
    src_dec_near = src_dec[Yary]
    src_flux_near = src_flux[Yary]
    src_name_near = src_name[Yary]
    src_map_near = np.vstack((src_flux_near,src_ra_near,src_dec_near))
    np.savetxt(outfilestr0+"integral.dat",src_map_sel.T,fmt="%s")
    return src_map_sel.T,ralb,rahb,declb,dechb,ra_centre,dec_centre,src_name_sel,src_num,src_map_near.T,src_name_near


if __name__ == "__main__":
    srcs,ralb,rahb,declb,dechb,ra_centre,dec_centre = src_map("config_srcmap.xml")
    print("srcs:::::::",srcs,ralb,rahb,declb,dechb)
    #srcs=srcs[srcs[:,0]>2.0]
    print(srcs,srcs.shape)
    mesrc =  np.loadtxt("./ME_LC/srclist.dat")
    mesrcerr =  np.loadtxt("./ME_LC/srclist_err.dat")
    lesrc =  np.loadtxt("./LE_LC/srclist.dat")
    lesrc_witherr =  np.loadtxt("./LE_LC/srclist_witherr.dat")
    lesrc_guan =  np.loadtxt("./LE_LC/srclist_guanju.dat")
    plt.scatter(lesrc[:,0],lesrc[:,1],marker='o',c = 'k',s=20*(lesrc[:,2]),edgecolor='k')
    plt.scatter(lesrc_witherr[:,0],lesrc_witherr[:,2],marker='o',c = 'r',s=12*(lesrc_witherr[:,4]),edgecolor='r')
    plt.scatter(lesrc_guan[:,0],lesrc_guan[:,1],marker='o',c = 'g',s=8*(lesrc_guan[:,2]),edgecolor='g')
    #plt.scatter(mesrcerr[:,0],mesrcerr[:,1],marker='o',c = 'g',s=50*(mesrcerr[:,2]))
    #plt.scatter(srcs[:,1],srcs[:,2],marker='o',s=30*(srcs[:,0]),c = 'b')
    plt.scatter(srcs[:,1],srcs[:,2],marker='o',s=20,c = 'b',edgecolor = 'b')
    att = np.loadtxt("att.dat")
    plt.plot(att[:,0],att[:,1],'r-')
    plt.grid()
    plt.show()


