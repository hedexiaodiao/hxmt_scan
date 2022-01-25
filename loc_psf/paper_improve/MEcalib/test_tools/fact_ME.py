#!/home/hxmt/nangyi/anaconda2/bin/python
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

#scnid = sys.argv[1]

fpath = '/hxmt/work/HXMT_scan_data/HE/nyRt2/'
foldertp = os.listdir(fpath)
#list1 = range(29500102,29500108,1)
#list0 = range(29500202,29500208,1)
#list2 = range(29500602,29500608,1)
#list3 = range(29500802,29500808,1)
#list4 = range(29500902,29500908,1)
#list5 = range(29501002,29501008,1)
#list6 = range(29501107,29501101,-1)
#list7 = range(29501207,29501201,-1)
#list3.remove(29500804)
#listall = list0 + list1 + list2 + list3 + list4 + list5 + list6 + list7

##########################################################
def gen_pol(scnid):
    alfa1=[]
    alfa2=[]
    alfa3=[]

    err1 = []
    err2=[]
    err3=[]

    norms1=[]
    norms2=[]
    norms3=[]

    factpath = "/hxmt/work/HXMT_scan_data/psfcrt/psf_factor/"
    fpath = "/hxmt/work/HXMT_scan_data/psfcrt/ME/"
    listp = os.listdir(fpath)
    scnid="ALL"
    #listp.remove(scnid)
    for i in os.listdir(fpath):
        if i[:5]!="P0101":
            continue
        dat1 = np.load(fpath+'%s/'%i + 'Angle_normT_0.npy')
        dat2 = np.load(fpath+'%s/'%i + 'Angle_normT_1.npy')
        dat3 = np.load(fpath+'%s/'%i + 'Angle_normT_2.npy')
        if len(alfa1)==0:
            frag = 6
            for j in range(frag):
                alfa1.append([itm for itm in dat1[0][j:len(dat1[0]):frag] if itm!=-99])
                norms1.append([itm for itm in dat1[1][j:len(dat1[0]):frag] if itm!=-99])
                err1.append([itm for itm in dat1[2][j:len(dat1[0]):frag] if itm!=-99])
            frag = 26
            for j in range(frag):
                alfa2.append([itm for itm in dat2[0][j:len(dat2[0]):frag] if itm!=-99])
                norms2.append([itm for itm in dat2[1][j:len(dat2[0]):frag] if itm!=-99])
                err2.append([itm for itm in dat2[2][j:len(dat2[0]):frag] if itm!=-99])
            frag = 6
            for j in range(frag):
                alfa3.append([itm for itm in dat3[0][j:len(dat3[0]):frag] if itm!=-99])
                norms3.append([itm for itm in dat3[1][j:len(dat3[0]):frag] if itm!=-99])
                err3.append([itm for itm in dat3[2][j:len(dat3[0]):frag] if itm!=-99])
        else:
            frag = 6
            for j in range(frag):
                alfa1[j].extend([itm for itm in dat1[0][j:len(dat1[0]):frag] if itm!=-99])
                norms1[j].extend([itm for itm in dat1[1][j:len(dat1[0]):frag] if itm!=-99])
                err1[j].extend([itm for itm in dat1[2][j:len(dat1[0]):frag] if itm!=-99])
            frag = 26
            for j in range(frag):
                alfa2[j].extend([itm for itm in dat2[0][j:len(dat2[0]):frag] if itm!=-99])
                norms2[j].extend([itm for itm in dat2[1][j:len(dat2[0]):frag] if itm!=-99])
                err2[j].extend([itm for itm in dat2[2][j:len(dat2[0]):frag] if itm!=-99])
            frag = 6
            for j in range(frag):
                alfa3[j].extend([itm for itm in dat3[0][j:len(dat3[0]):frag] if itm!=-99])
                norms3[j].extend([itm for itm in dat3[1][j:len(dat3[0]):frag] if itm!=-99])
                err3[j].extend([itm for itm in dat3[2][j:len(dat3[0]):frag] if itm!=-99])

    polt = [[],[],[]]
    for i in range(len(alfa3)):
        plt.errorbar(alfa3[i],np.array(norms3[i])/100.,yerr=np.array(err3[i])/100.,fmt='o',alpha=0.7)
        yesq = np.nonzero(np.array(norms3[i])>0)[0]
        pol = (np.polyfit(np.array(alfa3[i])[yesq], np.array(norms3[i])[yesq]/100., 7,w=100./np.array(err3[i])[yesq]))
        polt[0].append(pol)
        pol = np.poly1d(pol)
        plt.plot(np.arange(-5.7,5.7,0.05),pol(np.arange(-5.7,5.7,0.05)),'-')
        plt.ylim(0,2)
        plt.xlim(-5,5)
        plt.savefig(factpath+'%s_3ME%s.jpg'%(scnid,i))
        #plt.show()
        plt.close('all')

    for i in range(len(alfa1)):
        plt.errorbar(alfa1[i],np.array(norms1[i])/100.,yerr=np.array(err1[i])/100.,fmt='o',alpha=0.7)
        yesq = np.nonzero(np.array(norms1[i])>0)[0]
        pol = (np.polyfit(np.array(alfa1[i])[yesq], np.array(norms1[i])[yesq]/100., 7,w=100./np.array(err1[i])[yesq]))
        polt[2].append(pol)
        pol = np.poly1d(pol)
        plt.plot(np.arange(-5.7,5.7,0.05),pol(np.arange(-5.7,5.7,0.05)),'-')
        plt.ylim(0,2)
        plt.xlim(-5,5)
        plt.savefig(factpath+'%s_1ME%s.jpg'%(scnid,i))
        #plt.show()
        plt.close('all')

    for i in range(len(alfa2)):
        plt.errorbar(alfa2[i],np.array(norms2[i])/100.,yerr=np.array(err2[i])/100.,fmt='o',alpha=0.7)
        yesq = np.nonzero(np.array(norms2[i])>0)[0]
        pol = (np.polyfit(np.array(alfa2[i])[yesq], np.array(norms2[i])[yesq]/100.,7,w=100./np.array(err2[i])[yesq]))
        polt[1].append(pol)
        pol = np.poly1d(pol)
        plt.plot(np.arange(-5.7,5.7,0.05),pol(np.arange(-5.7,5.7,0.05)),'-')
        plt.ylim(0,2)
        plt.xlim(-2,2)
        plt.savefig(factpath+'%s_2ME%s.jpg'%(scnid,i))
        plt.close('all')

    np.save(factpath+'ME%s_polycrt.npy'%scnid,polt)

fpath = "/hxmt/work/HXMT_scan_data/psfcrt/ME/"
listp = os.listdir(fpath)
for i in os.listdir(fpath):
    if i[:5]!="P0101":
        continue
    gen_pol(i)

