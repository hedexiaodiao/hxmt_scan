#!/hxmt/soft/Develop/anaconda2/bin/python
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS,GCRS,ITRS
from astropy import coordinates as coord 
from astropy import units as u 
from astropy import time 
from astropy.time import Time
import astropy.io.fits as pf
import numpy as np
import sys
sys.path.append('/home/hxmt/lick/mylib/pylib')
from readxml import *

def fix_eci(infile,outfile):

   ###################read file#######################################     
    hdulist = pf.open(infile)
    orbit = hdulist[1].data
    ttime = orbit["Time"]
    mjd = (orbit.Time-3)/86400.+55927
    mjdtime = Time(mjd,format='mjd')
    utctime = mjdtime.isot
    itrs = ITRS(x=orbit["X"]*u.m, y=orbit["Y"]*u.m, z=orbit["Z"]*u.m,d_x=orbit["Vx"]*u.m/u.s,\
            d_y=orbit["Vy"]*u.m/u.s, d_z=orbit["Vz"]*u.m/u.s,representation=coord.CartesianRepresentation,\
            differential_cls=coord.CartesianDifferential,obstime=utctime)
    gcrs = itrs.transform_to(GCRS(obstime=utctime))
    gcrs.representation = 'cartesian'
    #######write to fits################################
    totcounts = gcrs.x.value.shape[0]
    phdu = pf.PrimaryHDU()
    atime = orbit.field(0)
    aorbit = np.zeros(totcounts)
    #a1 = gcrs.x.value
    #a2 = gcrs.y.value
    #a3 = gcrs.z.value
    #a4 = gcrs.v_x.value*1000
    #a5 = gcrs.v_y.value*1000
    #a6 = gcrs.v_z.value*1000
    a1 = gcrs.x.value*100
    a2 = gcrs.y.value*100
    a3 = gcrs.z.value*100
    a4 = gcrs.v_x.value*100000
    a5 = gcrs.v_y.value*100000
    a6 = gcrs.v_z.value*100000
    
    c1 = pf.Column(name="TIME",format="J",unit='s',array=atime)
    c2 = pf.Column(name="ORBIT",format="I",array=aorbit)
    #c3 = pf.Column(name="X",format="J",unit='m',bscale=0.01,array=a1)
    #c4 = pf.Column(name="Y",format="J",unit='m',bscale=0.01,array=a2)
    #c5 = pf.Column(name="Z",format="J",unit='m',bscale=0.01,array=a3)
    #c6 = pf.Column(name="Vx",format="J",unit='m/s',bscale=0.01,array=a4)
    #c7 = pf.Column(name="Vy",format="J",unit='m/s',bscale=0.01,array=a5)
    #c8 = pf.Column(name="Vz",format="J",unit='m/s',bscale=0.01,array=a6)
    #c9 = pf.Column(name="Lon",format="J",unit='degree',bscale=0.001,array=orbit["Lon"])
    #c10 = pf.Column(name="Lat",format="J",unit='degree',bscale=0.001,array=orbit["Lat"])
    #c11 = pf.Column(name="Alt",format="J",unit='m',bscale=0.001,array=orbit["Alt"])
    c3 = pf.Column(name="X",format="J",unit='m',array=a1)
    c4 = pf.Column(name="Y",format="J",unit='m',array=a2)
    c5 = pf.Column(name="Z",format="J",unit='m',array=a3)
    c6 = pf.Column(name="Vx",format="J",unit='m/s',array=a4)
    c7 = pf.Column(name="Vy",format="J",unit='m/s',array=a5)
    c8 = pf.Column(name="Vz",format="J",unit='m/s',array=a6)
    #c9 = pf.Column(name="Lon",format="J",unit='degree',array=orbit["Lon"])
    #c10 = pf.Column(name="Lat",format="J",unit='degree',array=orbit["Lat"])
    #c11 = pf.Column(name="Alt",format="J",unit='m',array=orbit["Alt"])
    c9 = pf.Column(name="Lon",format="J",unit='degree',array=1000*orbit["Lon"])
    c10 = pf.Column(name="Lat",format="J",unit='degree',array=1000*orbit["Lat"])
    c11 = pf.Column(name="Alt",format="J",unit='m',array=1000*orbit["Alt"])
    c12 = pf.Column(name="T",format="19A",array=aorbit)
    c13 = pf.Column(name="B",format="19A",array=aorbit)
    
    
    cols = pf.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    
    
    hdf = pf.open('/home/hxmt/lick/mylib/files/HXMT_P010129416301_Orbit_FFFFFF_V1_L1P.FITS')
    
    he_h0 = hdf[1].header
    #header0 = modify_header(he_h0,i,tstart,tstop)
    header0 = he_h0.copy()
    tbhdu.header.update(header0.cards[7:])
    phdu = pf.PrimaryHDU()
    tblist = []
    tblist.append(phdu)
    tblist.append(tbhdu)
    hdulist2 = pf.HDUList(tblist)
    hdulist2.writeto(outfile,overwrite=True)
if __name__ == "__main__":
    if  len(sys.argv)>2:
        infile = sys.argv[1]
        outfile = sys.argv[2]
        fix_eci(infile,outfile)
    else:
        print " Format:   python fix_eci_1N.py infile  outfile \n pls try again"
