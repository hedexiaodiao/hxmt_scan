#!/hxmt/home/lick/soft/anaconda2/bin/python
import sys
from astropy.io import fits as pf
import numpy as np

def loadsrclog(filename):
    norm = []
    num = []
    ra = []
    dec = []
    snr = []
    jd = 0
    error=[]
    name=[]
    chi = []
    pct=[]
    with open(filename,'r') as f:
        for line in f:
            #print jd
            if jd==0:
                if line.split()[0][0]=='P':
                    head = line.split()[1]
                    i = float(head)
                    if i not in num:
                        num.append(i)
                        #print i,num,jd,ra,line.split()[2],line.split()[3]
                        name.append((line.split('\t')[3][:-1])+"\t"+(line.split('\t')[2]))
                        jd = 1
            else:
                if jd > 0:
                    if jd>1:
                        if jd>2:
                            if jd>3:
                                if jd>4:
                                    #print line.split()
                                    pct.append(line.split()[2])
                                    jd=0
                                    #print pct
                                else:
                                    dec.append(float(line.split()[1]))
                                    jd=5
                                    #print dec
                            else:
                                ra.append(float(line.split()[1]))
                                jd =4
                                #print ra
                        else:
                            snr.append(float(line.split()[1]))
                            chi.append(float(line.split()[3]))
                            jd=3
                            #print chi
                    else:
                        norm.append(float(line.split()[1]))
                        #print norm
                        error.append(float(line.split()[5])/2 - float(line.split()[3])/2)
                        jd =2
            #raw_input("haha")
    #print 'haha'
    print(len(ra),len(snr))
    return ra,dec,norm,snr,error,name,num,chi,pct

def loadsrcnew(filename):
    norm = []
    num = []
    ra = []
    dec = []
    snr = []
    jd = 0
    error=[]
    raerr=[]
    decerr=[]
    near = []
    chi = []
    pct=[]
    with open(filename,'r') as f:
        for line in f:
            #print jd
            if jd==0:
                if line.split()[0][0]=='P':
                    head = line.split()[1]
                    #print head
                    i = float(head)
                    if i not in num:
                        num.append(i)
                        jd = 1
            else:
                #print "Good!" 
                if jd > 0:
                    if jd>1:
                        if jd>2:
                            if jd>3:
                                if jd>4:
                                    if jd>5:
                                        pct.append(line.split()[2])
                                        jd=0
                                    else:
                                        near.append(line)
                                        jd=6
                                else:
                                    dec.append(float(line.split()[1]))
                                    decerr.append(float(line.split()[5])/2 - float(line.split()[3])/2)
                                    jd=5
                            else:
                                ra.append(float(line.split()[1]))
                                raerr.append(float(line.split()[5])/2 - float(line.split()[3])/2)
                                jd =4
                        else:
                            snr.append(float(line.split()[1]))
                            chi.append(float(line.split()[3]))
                            jd=3
                    else:
                        norm.append(float(line.split()[1]))
                        error.append(float(line.split()[5])/2 - float(line.split()[3])/2)
                        jd=2
            #raw_input("haha")
    #print len(ra),len(near),"\n",near 
    return ra,raerr,dec,decerr,norm,error,snr,near,num,chi,pct

def findsrc(srcname):
    src_list = pf.open("/hxmt/work/HXMT_scan_data/HE/lick_tools/Srcs_gnl_swift.fits")
    dt = src_list[1].data.field(0)
    numid = np.nonzero(dt==srcname)[0][0]+1
    return numid
def gensrctxt(filepath,outpath):
    try:
        with open(filepath+"GAL_he_small.src",'r') as f:
            line = f.readline()
            fname = line.split()[0][5:]    
        ra,dec,norm,snr,error,name,num,chi,pct = loadsrclog(filepath+"GAL_he_small.src")
        with open(outpath+"src_P%s.txt"%fname,'w') as f:
            for i in range(len(ra)):
                ##print ('%s'%(snr[i]))
                f.write("%i\t"%(num[i]) + name[i] + "\t%s"%(ra[i]) + "\t%s"%(dec[i]) + "\t%s"%(norm[i]) + "\t%s"%(error[i]) + "\t%s"%(snr[i])+"\n")
    except:
        pass
    try:
        ra,raerr,dec,decerr,norm,error,snr,near,num,chi,pct = loadsrcnew(filepath+"GAL_he_smallNew.src")
        with open(outpath+"newsrclist.txt",'w') as f:
            for i in range(len(ra)):
                f.write("%i"%(num[i]) + "\t%s"%(ra[i]) + "\t%s"%(raerr[i]) + "\t%s"%(dec[i]) + "\t%s"%(decerr[i]) + "\t%s"%(norm[i]) + "\t%s"%(error[i]) + "\t%s"%(snr[i]) +"\t%s"%(chi[i])+"\t%s"%pct[i]+"\t%s"%(near[i]))
        for i in range(len(ra)):
            if snr[i]>3:
                if len(near[i].split('\t'))==6:
                    if float(near[i].split()[-1])<0.15:
                        numid = findsrc(near[i].split('\t')[1])
                        with open(outpath+"src_P%s.txt"%fname,'a') as f2:
                            f2.write("%i"%(num[i])+'\t'+str(numid)+'\t'+near[i].split('\t')[1]+'\t'+near[i].split('\t')[2]+'\t'+near[i].split('\t')[3]+'\t%s'%(norm[i]) + "\t%s"%(error[i]) + "\t%s"%(snr[i])+"\n")
    except:
        pass


if __name__ == "__main__":
    if len(sys.argv)<2:
        filename = raw_input("lack filename,please input:")
    gensrctxt(filename,filename)
