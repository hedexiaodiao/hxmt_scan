from astropy.io import fits as pf
import matplotlib.pyplot as plt
import numpy as np
def gti(time,peroid):
    time=np.array(time)
    if len(time)>0:
        t1 = time[:-1]
        t2 = time[1:]
        t_mark = t2-t1
        gtistop = np.append(t1[t_mark>peroid],time[-1])
        gtistart = np.append(time[0],t2[t_mark>peroid])
    else:
        gtistart=np.array([0]);gtistop=np.array([0])
    return gtistart,gtistop

#attfile='/hxmt/work/HXMT-DATA/1L/A01/P0101294/P0101294001/ACS/HXMT_P0101294001_Att_FFFFFF_V1_L1P.FITS'
def motionless(attfile):
    att=pf.open('%s'%attfile)
    time=att[1].data.field(0)
    ra=att[1].data.field('ra')
    dec=att[1].data.field('dec')
    #####take a point every 15 seconds
    time_peroid=time[-1]-time[0]
    point_num=time_peroid/15
    cut=int(len(time)/point_num)
    jumpt=time[::cut];jumpra=ra[::cut];jumpdec=dec[::cut]
    ########find motionless point
    jumpraa=abs(jumpra[1:]-jumpra[:-1])
    jumpdecc=abs(jumpdec[1:]-jumpdec[:-1])
    slecttime=[]
    for i in range(5,len(jumpraa)):
        if (len(jumpraa[i-5:i][jumpraa[i-5:i]<0.03])>=4)&(len(jumpraa[i-5:i][jumpdecc[i-5:i]<0.03])>=4):
            slecttime=np.append(slecttime,jumpt[i])    
    slectstart,slectstop=gti(slecttime,16)
    motionlesspoint=[];movepoint=[]
    if len(slecttime)>0:
        slectstart=slectstart-80
        slectstop=slectstop+60
        for j in range(len(slectstart)):
            motionlesspoint=np.append(motionlesspoint,time[(time>slectstart[j])&(time<slectstop[j])])
    movepoint=time[np.in1d(time,motionlesspoint)==0]
    gtistart,gtistop=gti(movepoint,0.5)

    return gtistart,gtistop


def find_slew(attfile0):
   atttable=pf.open(attfile0)[1].data
   atttm=atttable.field(0)
   attra=atttable.field('ra');attdec=atttable.field('dec')
   hhra=attra.copy();hhdec=attdec.copy()
   hhra[hhra>180]=hhra[hhra>180]-360
   select_time1=[];select_time2=[];
   tstart,tstop=gti(atttm,4)
   print len(tstart)
   for num in range(len(tstart)):
     mask=(atttm>=tstart[num])&(atttm<=tstop[num])
     time=atttm[mask]
     interval1=np.sqrt((attra[mask][1:]-attra[mask][:-1])**2+(attdec[mask][1:]-attdec[mask][:-1])**2)
     speed1=interval1/(time[1:]-time[:-1])
     interval2=np.sqrt((hhra[mask][1:]-hhra[mask][:-1])**2+(hhdec[mask][1:]-hhdec[mask][:-1])**2)
     speed2=interval2/(time[1:]-time[:-1])
     mask1=(speed1>0.023)&(speed1<0.13)
     mask2=(speed2>0.023)&(speed2<0.13)
     select_time1=np.append(select_time1,time[1:][mask1])
     select_time2=np.append(select_time2,time[1:][mask2])
#     print num,len(atttm),len(select_time1),len(select_time2)
     if len(select_time1)>len(select_time2):
        select_time=select_time1
     else:
        select_time=select_time2
     print num,len(atttm),len(select_time1),len(select_time2)
   return select_time,attra[np.in1d(atttm,select_time)],attdec[np.in1d(atttm,select_time)]




