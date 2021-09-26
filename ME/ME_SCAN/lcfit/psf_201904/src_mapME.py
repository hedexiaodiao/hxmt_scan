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
import sys
#sys.Wpath.append("/hxmt/home/saina/tar/tools/")
import astropy.io.fits as pf
import numpy as np
from math import *
from readxml import *
import time
import Quat as quat
#import matplotlib.pyplot as plt
#############read ra dec for scan stripe################
def src_map(cfg,rac,decc,slist):
    readcfg = loadDom(cfg)
    infilestr = readcfg.getTagText("infilenamelist")
    inpathstr = readcfg.getTagText("inpath")
    outpathstr = readcfg.getTagText("outpath")
    inpathstr = inpathstr.strip()
    outpathstr = outpathstr.strip()
    evtfilestr = infilestr.split()[0]
    evtfilestr2 = infilestr.split()[-2]
    infilestr2 = infilestr.split()[-1]
    infile = (infilestr2).encode()###
    evtfile = (evtfilestr).encode()###
    evtfile2 = (evtfilestr2).encode()###
    outpath=outpathstr.encode()
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
    etime = etime1
    print "tstart,tstop: ",tstart,tstop
    tb3sel = tb3[(tb3.field(0)<=tstop) & (tb3.field(0)>=tstart)]
    print tb3sel.shape
    time.sleep(5)
    q1_list = tb3.field(1)[np.in1d(qtime,etime)]
    q2_list = tb3.field(2)[np.in1d(qtime,etime)]
    q3_list = tb3.field(3)[np.in1d(qtime,etime)]
    stripe_ra_list = []
    stripe_dec_list = []
    print "q list size: ",q1_list.shape
    fatt=open(outpath+"att.dat",'w')
    fq=open(outpath+"q.dat",'w')
    for i in range(0,q1_list.shape[0]):
        qu = [q1_list[i],q2_list[i],q3_list[i]]
        qu = qu+[sqrt(fabs(1-qu[0]*qu[0]-qu[1]*qu[1]-qu[2]*qu[2]))]
        quat1 = quat.Quat(quat.normalize(qu))
        stripe_ra_list.append(quat1.ra)
        stripe_dec_list.append(quat1.dec)
    time.sleep(5)
    fatt.close()
    fq.close()
    stripe_ra_list = np.array(stripe_ra_list)
    stripe_dec_list = np.array(stripe_dec_list)
    ra_centre = stripe_ra_list.mean()       
    dec_centre = stripe_dec_list.mean()
    ralb = stripe_ra_list.min()
    rahb = stripe_ra_list.max()
    declb = stripe_dec_list.min()
    dechb = stripe_dec_list.max()
    #src_list = pf.open("/hxmt/work/HDPC/KS_Test/xspec_sim_hxmt_survey/s_map/gnrl_refr_cat_0031.fits")
    #src_map = src_list[1].data
    #src_name = src_map.field(2)
    #src_ra = src_map.field(4)
    #src_dec = src_map.field(5)
    #src_flux = src_map.field(28)
    #src_list = pf.open("/hxmt/work/HXMT_scan_data/psfcrt/tools/Srcs_IGR_SWIFT.fits")
    src_list = pf.open(slist)
    src_map = src_list[1].data
    src_name = src_map.field(0)
    src_ra = src_map.field(1)
    src_dec = src_map.field(2)
    src_flux = src_map.field(4)
    Yary=[]
    num_list=[]
    for i,itra in enumerate(src_ra):
        if ((((itra - stripe_ra_list)**2 + (src_dec[i]-stripe_dec_list)**2)<=17)).sum()>0:
            Yary.append(1)
            num_list.append(i)
        else:
            Yary.append(0)
    print sum(Yary)
    Yary=(np.array(Yary)>0)
    src_num = np.nonzero(Yary)[0]
    #radius2r = 20**2
    #radius2d = 14**2
    #src_ra_sel = src_ra[Yary & ((src_ra[:]-rac)**2<=radius2r) & ((src_dec[:]-decc)**2<=radius2d)]
    #src_dec_sel = src_dec[Yary & ((src_ra[:]-rac)**2<=radius2r) & ((src_dec[:]-decc)**2<=radius2d)]
    #src_flux_sel = src_flux[Yary & ((src_ra[:]-rac)**2<=radius2r) & ((src_dec[:]-decc)**2<=radius2d)]
    #src_name_sel = src_name[Yary & ((src_ra[:]-rac)**2<=radius2r) & ((src_dec[:]-decc)**2<=radius2d)]
    #src_num_sel = num_list
    src_ra_sel = src_ra[Yary]
    src_dec_sel = src_dec[Yary]
    src_flux_sel = src_flux[Yary]
    src_name_sel = src_name[Yary]
    src_flux_tot = src_flux_sel
    #print src_ra_sel,src_dec_sel,src_name_sel
    print len(src_ra_sel),len(src_flux_sel)
    src_map_sel_name = np.vstack((src_flux_tot,src_ra_sel,src_dec_sel,src_name_sel))
    src_map_sel = np.vstack((src_flux_tot,src_ra_sel,src_dec_sel))
    print src_map_sel.shape,len(src_ra_sel),len(src_dec_sel),len(src_flux_sel),len(src_name_sel)
    np.savetxt(outpath+"maxi&integral&swift.dat",src_map_sel.T,fmt="%s")
    np.savetxt(outpath+"maxi&integral&swift_name.dat",src_map_sel_name.T,fmt="%s")
    return src_map_sel.T,ralb,rahb,declb,dechb,ra_centre,dec_centre,src_name_sel

'''if __name__ == "__main__":
    srcs,ralb,rahb,declb,dechb,ra_centre,dec_centre,src_name,src_num = src_map("config_me_small.xml")
    #print "srcs:::::::",srcs,ralb,rahb,declb,dechb
    #print srcs,srcs.shape
    #plt.scatter(srcs[:,1],srcs[:,2],marker='o',s=10,c = 'b',edgecolor = 'b')
    att = pf.open("./att.fits")
    ra=att[1].data.field(1)
    dec=att[1].data.field(2)
    #plt.plot(ra,dec,'r-')
    #plt.grid()
    #plt.show()
'''
    
