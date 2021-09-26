#!/hxmt/home/lick/soft/anaconda2/bin/python
'''
Copyright 2016, HXMT Ground System
All right reserved
File Name: loadmod.py for sample

.......

Author: lick
Version: 0.1
Date:  2016-12-11
Email: lick@ihep.ac.cn
......
'''
from astropy.io import fits as pf
import sys
import os
import time
import numpy as np
import math

from exec_cmd import *

def cpha(crf,fake):
    fname = crf.split('/')[-1][9:13]
    print('crf,fake,fname:',crf,fake,fname)
    hd = pf.open(crf)
    tb = hd[1].data
    t = tb.field(0)
    t0  = t[0]
    cr2 = tb.field(1)
    err2 = tb.field(2)
    ran = time.time()
    f = open("%stemppha%s.txt"%(fname,ran),"w")
    f2 = open("%stempbad%s.dat"%(fname,ran),"w")
    for i,r in enumerate(cr2):
        tst = i 
        tsp = i+1
        r = r
        if i ==0:
            if r <1:
                tl = str(i)+" "
                f2.writelines(tl)
        if i >0 and i <len(cr2)-1:
            if r<1:
                if cr2[i-1]>1:
                    tl2 = str(i)+" "
                    f2.writelines(tl2)
            if r>1:
                if cr2[i-1]<1:
                    tl3 = str(i-1)+"\n"
                    f2.writelines(tl3)
        if i == len(cr2)-1:
            if r<1:
                if cr2[i-1]<1:
                    tl4 =str(i)+"\n"
                    f2.writelines(tl4)
                if cr2[i-1]>=1:
                    tl4 =str(i)+" "+str(i)+"\n"
                    f2.writelines(tl4)
        re = err2[i]
        line = str(tst)+" "+str(tsp)+" "+str(r)+" "+str(re)+"\n"
        f.writelines(line)
    for i,r in enumerate(err2):
        tst = i 
        tsp = i+1
        r = r
        if i ==0:
            if r <0.1:
                tl = str(i)+" "
                f2.writelines(tl)
        if i >0 and i <len(err2)-1:
            if r<0.1:
                if err2[i-1]>0.1:
                    tl2 = str(i)+" "
                    f2.writelines(tl2)
            if r>0.1:
                if err2[i-1]<0.1:
                    tl3 = str(i-1)+"\n"
                    f2.writelines(tl3)
        if i == len(err2)-1:
            if r<0.1:
                if err2[i-1]<0.1:
                    tl4 =str(i)+"\n"
                    f2.writelines(tl4)
                if err2[i-1]>=0.1:
                    tl4 =str(i)+" "+str(i)+"\n"
                    f2.writelines(tl4)
        #print line
    f.close()
    f2.close()
    my_env = os.environ
    commond = 'flx2xsp %stemppha%s.txt %s.pha %s.rsp' % (fname, ran, fake, fake)
    cmd_tem = 'export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;' + commond
    exec_code, exec_log = exec_cmd(cmd_tem, my_env)
    print(cmd_tem)
    print('exec_code:', exec_code)

    os.system('rm %stemppha%s.txt'%(fname,ran))

    commond = 'grppha %s.pha !%s.grp "bad %stempbad%s.dat&exit"'%(fake,fake,fname,ran)
    cmd_tem = 'export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;' + commond
    exec_code, exec_log = exec_cmd(cmd_tem, my_env)
    print(cmd_tem)
    print('exec_code:', exec_code)
    #f3=open('%stempbad%s.dat'%(fname,ran))
    #print f3.readlines()
    os.system('rm %stempbad%s.dat'%(fname,ran))
    print("generation of the faked pha file and rsp file : success")

if __name__ == "__main__":
    print("test")
    cpha("/hxmt/work/USERS/liuyuan/wangchen/tar/P1294/p0101294510/norm/GAL_LE_small0.pha","fake")
