#from multiprocessing import Pool
import sys,os
import numpy as np
print('success3.00')
from readxml import *
print('success3.01')
from xspec import *
import numpy as np
from astropy.io import fits as pf
print('success3.02')
import funchxmt_paraodic as hxmtpsf
print('success3.03')
import time

'''
def line(a,b,dist):
    x_mid = (a[0]+b[0])/2.0
    y_mid = (a[1]+b[1])/2.0
    k_line = (a[1] - b[1]) / (a[0] - b[0]) if a[0] != b[0] else None
    #print k_line,x_mid,y_mid
    X,Y = symbols("X,Y")
    if k_line == 0:
        relt = [[x_mid,y_mid-dist],[x_mid,y_mid+dist]]
    elif k_line == None:
        relt = [[x_mid-dist,y_mid],[x_mid+dist,y_mid]]
    else:
        k_line = -1.0/k_line
        relt = solve((k_line * X + y_mid - k_line * x_mid - Y,X ** 2 - 2 * X * x_mid + x_mid ** 2 + Y ** 2 - 2 * Y * y_mid + y_mid ** 2 - dist**2),X,Y)
    return np.array(relt)
'''
def Sel_Median(array,sigma = 2.5):
    x = np.copy(array)
    od = np.arange(len(x))
    while 1:
        md = np.median(x)
        td = np.std(x)
        bad = x > md + sigma*td
        if bad.sum()>0:
            x = x[~bad]
            od = od[~bad]
        else:
            return x,od,x[0].mean()

def FindGap(inptm,gap=1):
    times = np.copy(inptm)
    a=((times[1:]-times[:-1])>gap)
    c=np.nonzero(a)[0]
    potstart = np.r_[0,c+1]
    potstop = np.r_[c,times.shape[0]-1]
    return potstart,potstop

def distance(ratpp,dectpp,rapnt,decpnt):
    ytp = (np.sin(decpnt*np.pi/360.0 - dectpp*np.pi/360.0))**2 + np.cos(decpnt*np.pi/180.0)*np.cos(dectpp*np.pi/180.0)*(np.sin(ratpp*np.pi/360.0 - rapnt*np.pi/360.0))**2
    ytp = np.arccos(1-2*ytp)*180./np.pi
    return ytp


def GetAttInf(attfile,lcfile):

    atthd = pf.open(attfile)
    atttb = atthd[1].data
    atttm = atttb.field(0)
    #rapnt = float(atthd[0].header['RA_OBJ'])
    #decpnt= float(atthd[0].header['DEC_OBJ'])
    #dtobs = atthd[0].header['DATE-OBS']
    lchd = pf.open(lcfile)
    lctb = lchd[1].data
    lctm = lctb.field(0)
    a=((lctm[1:]-lctm[:-1])!=1)
    c=np.nonzero(a)[0]
    pot=np.r_[0,c,lctm.shape[0]-1]
    c=np.zeros(len(atttm))
    c=(c==1)
    for i in range(len(pot)-1):
        c=(c|((atttm>=lctm[pot[i]+1])&(atttm<=lctm[pot[i+1]])))

    tbatt = atttb[c]
    ratp = tbatt.field(1)
    dectp = tbatt.field(2)

    return 1,1,1,ratp,dectp#rapnt,decpnt,dtobs,ratp,dectp

def GetXmlInf(xmlfile):
    readcfg = loadDom(xmlfile)
    infilestr = readcfg.getTagText("infilenamelist")
    inpathstr = readcfg.getTagText("inpath")
    outpathstr = readcfg.getTagText("outpath")
    outfilestr = readcfg.getTagText("outfilename")
    print(inpathstr,infilestr,outpathstr)
    inpathstr = inpathstr.strip()
    infilestr0 = infilestr.split()[0]
    infilestr1 = infilestr.split()[1]
    infilestr2 = infilestr.split()[2]
    infilestr3 = infilestr.split()[3]
    infile0 = (inpathstr+infilestr0)#.encode()
    infile1 = (inpathstr+infilestr1)#.encode()
    infile2 = (inpathstr+infilestr2)#.encode()
    infile3 = (inpathstr+infilestr3)#.encode()
    outpathstr = outpathstr.strip()
    outfilestr0 = outfilestr.split()[0]
    outfilePri = (outpathstr+outfilestr0)#.encode()
    infilePri = infile3.split("/")[-1][:13]
    return infile0,infile1,infile2,infile3,infilePri,outfilePri

from xspec import *
def InitModel_T(src,outfilePri):
    print("New pot Runing >>>>>>>>>>> %s >>>>>>>>>>>\n"%src[0])
    time1= time.time()
    #from xspec import *
    AllModels.addPyMod(hxmtpsf.hxmtpsf, hxmtpsf.myModelParInfo, 'add')
    AllModels.clear()
    AllData.clear()
    AllData("1:1 %s0.grp 2:2 %s1.grp 3:3 %s2.grp"%(outfilePri,outfilePri,outfilePri))
    d1 = AllData(1)
    d2 = AllData(2)
    d3 = AllData(3)
    AllModels += "powerlaw+hxmtpsf"
    m1 = AllModels(1)
    m2 = AllModels(2)
    m3 = AllModels(3)
    m1.setPars({1:"0 -1",2:"1. -0.1"})
    m2.setPars({1:"0 -1",2:"1. -0.1"})
    m3.setPars({1:"0 -1",2:"1. -0.1"})
    Fit.query='yes'
    ra,dec,norm = src[1],src[2],src[0]
    print(src,norm)
    m1.setPars({3:" 1          -1   "})
    m1.setPars({4:"60    -0.001"})
    m1.setPars({5:"%f,-0.01"%(ra)})
    m1.setPars({6:"%f,-0.01"%(dec)})
    m1.setPars({7:"%f,0.01"%(int(norm*18./100.) if int(norm*18./100.)<1000. else 900.)})
    m2.setPars({4:"0 -0.001"})
    m3.setPars({4:"-60  -0.001"})
    Fit.perform()
    Fit.error("max 100 1.0 1-7")
    tpchi = Fit.statistic
    tpdof = Fit.dof
    try:
        snr = 2*m1(7).values[0]/abs(m1(7).error[1]-m1(7).error[0])
        norm = m1(7).values[0]
    except:
        snr = -1
        norm = -1
    m1.setPars({7:"0.,-0.01"})
    zerochi = Fit.statistic
    deltachi = abs(tpchi - zerochi)
    low = m1(7).error[0]
    upp = m1(7).error[1]
    if snr!=-1:
        m1.setPars({7:"%f,0.01"%norm})
        m1(5).frozen=False
        m1(6).frozen=False
        Fit.perform()
        Fit.error("max 100 1.0 1-7")
        nowchi = Fit.statistic
        nowdof = Fit.dof
        print(nowchi,nowdof,tpchi,tpdof)
        f = Fit.ftest(nowchi,nowdof,tpchi,tpdof)
    else:
        f=-1
    AllModels.clear()
    AllData.clear()
    print("New pot END!!!!!!!!!!!! %s >>>>>>>>>>>\n"%src[0])
    time2 = time.time()
    del AllChains,AllData,AllModels,Fit,Plot,Xset
    return [ra,dec,norm,tpchi,tpdof,low,upp,f]

from xspec import *
def InitModel(src,outfilePri):
    print("New pot Runing >>>>>>>>>>> %s >>>>>>>>>>>\n"%src[0])
    time1= time.time()
    #from xspec import *
    #Model('grbm')
    AllModels
    #print(AllModels(1)(4))
    AllData.clear()
    AllModels.clear()
    AllModels.addPyMod(hxmtpsf.hxmtpsf, hxmtpsf.myModelParInfo, 'add')
    AllData("1:1 %s0.grp 2:2 %s1.grp 3:3 %s2.grp"%(outfilePri,outfilePri,outfilePri))
    d1 = AllData(1)
    d2 = AllData(2)
    d3 = AllData(3)
    AllModels += "powerlaw+hxmtpsf"
    m1 = AllModels(1)
    m2 = AllModels(2)
    m3 = AllModels(3)
    m1.setPars({1:"0 -1",2:"0. 0.1"})
    m2.setPars({1:"0 -1",2:"0. 0.1"})
    m3.setPars({1:"0 -1",2:"0. 0.1"})
    Fit.query='yes'
    ra,dec,norm = src[1],src[2],src[0]
    print(src,norm)
    m1.setPars({3:" 1          -1   "})
    m1.setPars({4:"60    -0.001"})
    m1.setPars({5:"%f,-0.01"%(ra)})
    m1.setPars({6:"%f,-0.01"%(dec)})
    m1.setPars({7:"%f,0.01"%(int(norm*18./100.) if int(norm*18./100.)<1000. else 900.)})
    m2.setPars({4:"0 -0.001"})
    m3.setPars({4:"-60  -0.001"})
    Fit.perform()
    Fit.error("max 100 1.0 1-7")
    tpchi = Fit.statistic
    tpdof = Fit.dof
    try:
        snr = 2*m1(7).values[0]/abs(m1(7).error[1]-m1(7).error[0])
        norm = m1(7).values[0]
    except:
        snr = -1
        norm = -1
    m1.setPars({7:"0.,-0.01"})
    zerochi = Fit.statistic
    deltachi = abs(tpchi - zerochi)
    low = m1(7).error[0]
    upp = m1(7).error[1]
    AllModels.clear()
    AllData.clear()
    print("New pot END!!!!!!!!!!!! %s >>>>>>>>>>>\n"%src[0])
    time2 = time.time()
    #del AllChains,AllData,AllModels,Fit,Plot,Xset
    return [ra,dec,norm,tpchi,tpdof,low,upp,deltachi]

def Sel_srcs(src):
    peaks = 0
    for i in [60,0,-60]:
        psf = hxmtpsf.Cal_psf([0,i,src[1],src[2]])
        a,b = FindGap(np.array(range(len(psf)))[psf>0])
        if (psf.sum()>0):
            peaks+=1
        if (len(a)>2):
            peaks+=1
    return peaks

def InitModel_fixbkg(src,outfilePri):
    print("New pot Runing >>>>>>>>>>> %s >>>>>>>>>>>\n"%src[0])
    time1= time.time()
    #from xspec import *
    AllModels.addPyMod(hxmtpsf.hxmtpsf, hxmtpsf.myModelParInfo, 'add')
    AllModels.clear()
    AllData.clear()
    AllData("1:1 %s0.grp 2:2 %s1.grp 3:3 %s2.grp"%(outfilePri,outfilePri,outfilePri))
    d1 = AllData(1)
    d2 = AllData(2)
    d3 = AllData(3)
    AllModels += "powerlaw+hxmtpsf"
    m1 = AllModels(1)
    m2 = AllModels(2)
    m3 = AllModels(3)
    m1.setPars({1:"0 -1",2:"%s -0.1"%(Sel_Median(d1.values)[2]-0.8)})
    m2.setPars({1:"0 -1",2:"%s -0.1"%(Sel_Median(d2.values)[2]-0.8)})
    m3.setPars({1:"0 -1",2:"%s -0.1"%(Sel_Median(d3.values)[2]-0.8)})
    Fit.query='yes'
    ra,dec,norm = src[1],src[2],src[0]
    print(src,norm)
    m1.setPars({3:" 1          -1   "})
    m1.setPars({4:"60    -0.001"})
    m1.setPars({5:"%f,-0.01"%(ra)})
    m1.setPars({6:"%f,-0.01"%(dec)})
    m1.setPars({7:"%f,0.01"%(int(norm*18./100.) if int(norm*18./100.)<1000. else 900.)})
    m2.setPars({4:"0 -0.001"})
    m3.setPars({4:"-60  -0.001"})
    Fit.perform()
    Fit.error("max 100 1.0 1-7")
    tpchi = Fit.statistic
    tpdof = Fit.dof
    try:
        snr = 2*m1(7).values[0]/abs(m1(7).error[1]-m1(7).error[0])
        norm = m1(7).values[0]
    except:
        snr = -1
        norm = -1
    m1.setPars({7:"0.,-0.01"})
    zerochi = Fit.statistic
    deltachi = abs(tpchi - zerochi)
    low = m1(7).error[0]
    upp = m1(7).error[1]
    AllModels.clear()
    AllData.clear()
    print("New pot END!!!!!!!!!!!! %s >>>>>>>>>>>\n"%src[0])
    time2 = time.time()
    del AllChains,AllData,AllModels,Fit,Plot,Xset
    return [ra,dec,norm,tpchi,tpdof,low,upp,deltachi]


def AddNewPot(ratp,dectp):
    ratpp= np.array([])
    dectpp= np.array([])
    stepadd = 150

    for index,ra in enumerate(ratp[::stepadd]):
        dec = dectp[::stepadd][index]
        try:
            ytp = distance(ratpp,dectpp,ra,dec)
            ytp = (ytp < 1.5).sum()
            if ytp <1:
                ratpp = np.append(ratpp,ra)
                dectpp = np.append(dectpp,dec)
        except:
            pass

    addds = [[2.5,0.],[0.,2.5],[-2.5,0.],[0.,-2.5],[3.6,3.6],[3.6,-3.6],[-3.6,3.6],[-3.6,-3.6]]
    for index,ra in enumerate(ratpp):
        dec = dectpp[index]
        for addpot in addds:
            addra = ra + addpot[0]
            if addra>360.:
                addra = addra-360.
            elif addra<0.:
                addra = addra+360.
            adddec = dec + addpot[1]
            ytp = distance(ratpp,dectpp,addra,adddec)
            ytp = (ytp < 3).sum()
            if ytp <1:
                ratpp = np.append(ratpp,addra)
                dectpp = np.append(dectpp,adddec)

    return ratpp,dectpp

def AddNewPot2(ratp,dectp):
    ratpp= np.array([])
    dectpp= np.array([])
    stepadd = 10

    for index,ra in enumerate(ratp[::stepadd]):
        dec = dectp[::stepadd][index]
        try:
            ytp = distance(ratpp,dectpp,ra,dec)
            ytp = (ytp < 2).sum()
            if ytp <1:
                ratpp = np.append(ratpp,ra)
                dectpp = np.append(dectpp,dec)
        except:
            pass
    addds = zip(sum(np.mgrid[-5.7:5.7:115j,-5.7:5.7:115j][0].tolist(),[]),sum(np.mgrid[-5.7:5.7:115j,-5.7:5.7:115j][1].tolist(),[]))
    #addds = [[2.5,0.],[0.,2.5],[-2.5,0.],[0.,-2.5],[3.6,3.6],[3.6,-3.6],[-3.6,3.6],[-3.6,-3.6]]
    for index,ra in enumerate(ratpp):
        dec = dectpp[index]
        for addpot in addds:
            addra = ra + round(addpot[0],2)
            if addra>360.:
                addra = addra-360.
            elif addra<0.:
                addra = addra+360.
            adddec = dec + round(addpot[1],2)
            ytp = distance(ratpp,dectpp,addra,adddec)
            ytp = (ytp < 0.1).sum()
            if ytp <1:
                ratpp = np.append(ratpp,addra)
                dectpp = np.append(dectpp,adddec)

    return ratpp,dectpp

def exit_opin(ramax,decmax,ramin,decmin,ratp,dectp):
    yep = (ratp > ramin-5.7)&(ratp < ramax+5.7)&(dectp > decmin-5.7)&(dectp < decmax+5.7)
    return np.sum(yep)


def divid_field(ramax,decmax,ramin,decmin):
    xdist = distance(ramin,decmin,ramax,decmin)
    ydist = distance(ramin,decmin,ramin,decmax)
    a=[ramin,ramax] if xdist<1 else [ramin,(ramin+ramax)*0.5,ramax]
    b=[decmin,decmax] if ydist<1 else [decmin,(decmin+decmax)*0.5,decmax]
    rt = []
    for i in xrange(len(a)-1):
        for j in xrange(len(b)-1):
            rt.append([round(a[i+1],1),round(b[j+1],1),round(a[i],1),round(b[j],1)])
    return rt


def priv_pot(ramax,decmax,ramin,decmin):
    a=np.arange(round(ramin,1),round(ramax,1)-0.01,0.2).tolist()
    b=np.arange(round(decmin,1),round(decmax,1)-0.01,0.2).tolist()
    addds = zip(sum(np.around(np.meshgrid(a,b),decimals=2)[0].tolist(),[]),sum(np.around(np.meshgrid(a,b),decimals=2)[1].tolist(),[]))
    return addds

def AddNewPot3(ratp,dectp):
    newpot = []
    ramax = max(ratp)
    ramin = min(ratp)
    if ramax - ramin > 180.:
        ratp[ratp<180] = ratp[ratp<180] + 360.
        ramax = max(ratp)
        ramin = min(ratp)
    decmax = np.median(dectp) + 30.
    decmin = np.median(dectp) - 30.
    ramax = np.median(ratp) + 30.
    ramin = np.median(ratp) - 30.
    yeps = distance(ramax,decmax,ramin,decmin) > 1.42
    fields = []
    fields += divid_field(ramax,decmax,ramin,decmin)
    fields_yes = [[ramax,decmax,ramin,decmin]]
    while yeps:
        fields_yes = []
        for field in fields:
            yep = exit_opin(field[0],field[1],field[2],field[3],ratp,dectp) > 0
            if yep:
                fields_yes += [field]
        #print ":::1",yep
        fields = []
        for field in fields_yes:
            fields += divid_field(field[0],field[1],field[2],field[3])
        #print ":::2",fields
        yeps = distance(field[0],field[1],field[2],field[3]) > 1.42
        #print ":::3",yeps,distance(field[0],field[1],field[2],field[3]),field[0],field[1],field[2],field[3]
        #print divid_field(field[0],field[1],field[2],field[3])
        #print '*****'
        #raw_input('::')

    for field in fields:
        newpot += priv_pot(field[0],field[1],field[2],field[3])
    #print fields[:5]
    newpot = np.array(newpot)
    return newpot[:,0],newpot[:,1]

### Old Version ###
''' 
def divid_field(ramax,decmax,ramin,decmin):
    la = int((ramax-ramin)*5/2.)
    lb = int((decmax-decmin)*5/2.)
    a=[ramin,ramin+0.2] if la==0 else [ramin,ramin+la*0.2,ramin+la*0.4]
    b=[decmin,decmin+0.2] if lb==0 else [decmin,decmin+lb*0.2,decmin+lb*0.4]
    adds = zip(sum(np.meshgrid(a,b)[0].tolist(),[]),sum(np.meshgrid(a,b)[1].tolist(),[]))
    rt = []
    la = 1 if la==0 else la
    lb = 1 if lb==0 else lb
    for i in adds:
        if (i[0]+la*0.2,i[1]+lb*0.2) in adds:
            rt.append([round(i[0]+la*0.2,1),round(i[1]+lb*0.2,1),round(i[0],1),round(i[1],1)])
    return rt
#    return [[round(ramax,1),round(decmax,1),(round(ramax/2.+ramin/2.,1)),round(decmax/2.+decmin/2.,1)],[round(ramax/2.+ramin/2.,1),round(decmax,1),round(ramin,1),round(decmax/2.+decmin/2.,1)],[round(ramax/2.+ramin/2.,1),round(decmax/2.+decmin/2.,1),round(ramin,1),round(decmin,1)],[round(ramax,1),round(decmax/2.+decmin/2.,1),round(ramax/2.+ramin/2.,1),round(decmin,1)]]
'''
