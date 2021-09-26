#!/hxmt/soft/Develop/anaconda2/bin/python
'''
......method of application.............
./loadmod.py config.xml ALL
........................................
'''
###############CLass and Funciton defined by Others##########################
from xspec import *
import math
import numpy as np
import time
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
import astropy.io.fits as pf
import matplotlib.gridspec as gridspec
time1= time.time()
Program_dir = '/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN'
print "fit begin: ", time1
if len(sys.argv)<2:
	print "Need the config file!"
	sys.exit(1)
print sys.argv[0],sys.argv[1]
cfg = sys.argv[1]
#cfg = 'config_me_small.xml'
sys.path.append("%s/lcfit/psf_201904/"%Program_dir)
###sys.path.append("/sharefs/hbkg/user/saina/tar/psf_201904/")
import src_mapME as map
from readxml import *
import createpha_betax as createpha
import hxmt_paraodic_test as hxmtpsf
print "import class success"
#os.system("mkdir result/")
#######read input output template paths and files from  config file##########
readcfg = loadDom(cfg)
infilestr = readcfg.getTagText("infilenamelist")
inpathstr = readcfg.getTagText("inpath")
outpathstr = readcfg.getTagText("outpath")
outfilestr = readcfg.getTagText("outfilename")
print inpathstr,infilestr,outpathstr
time.sleep(3)
inpathstr = inpathstr.strip()
inpath=inpathstr.encode()
infilestr0 = infilestr.split()[0]
infilestr1 = infilestr.split()[1]
infilestr2 = infilestr.split()[2]
infilestr3 = infilestr.split()[-1]
infile0 = (infilestr0).encode()###
infile1 = (infilestr1).encode()###(inpathstr+infilestr1).encode()
infile2 = (infilestr2).encode()###
infile3 = (infilestr3).encode()###
#filename = infile3.split("/")[-1][5:18]
outpathstr = outpathstr.strip()
outfilestr0 = outfilestr.split()[0]
outpath=outpathstr.encode()
outfile0 = (outpathstr+outfilestr0+".log").encode()###outpathstr+outfilestr0
outfile1 = (outpathstr+outfilestr0+".par").encode()
outfile2 = (outpathstr+outfilestr0+".eps").encode()
outfile3 = (outpathstr+outfilestr0).encode()
outfile4 = (outpathstr+outfilestr0+".src").encode()
outfile5 = (outpathstr+outfilestr0).encode()
###outfile5 = (outfilestr0.split("/")[-1]).encode()###inpathstr[:-3]+
print ("the count rate file : ",infile0,infile1,infile2,infile3)
print ("the out put file : ",outfile0,outfile1,outfile2)
time.sleep(5)
slist='%s/lcfit/psf_201904/Srcs_IGR_SWIFT_06.fits'%Program_dir
src_list = pf.open(slist)

src_map = src_list[1].data
src_name = src_map.field(0)
src_ra = src_map.field(1)
src_dec = src_map.field(2)
src_flux = src_map.field(4)

##############load the hxmtpsf model#################################
try:
	logFile = Xset.openLog(outfile0)
	parFile = open(outfile1,"w")
	srcFile = open(outfile4,"w") 
	newsrcFile = open(outfile4[:-4]+"New.src","w")
except:
	print "Open log or par file failed!"
	sys.exit(1)

print infile0,infile1,infile2,hxmtpsf.infile
AllModels.addPyMod(hxmtpsf.hxmtpsf, hxmtpsf.myModelParInfo, 'add')

print "add hxmtpsf to the model"
time.sleep(5)
#########calculate the src map in the region of interesting ################
crab_unit = 1000
bkg_flux = 25
flux_thres = 0.0
centre_me = np.loadtxt(outpath+"centre.dat")
rac=centre_me[1]
decc=centre_me[2]
srcs,ralb,rahb,declb,dechb,ra_centre,dec_centre,src_name_all = map.src_map(cfg,rac,decc,slist)
print "srcs:::::::",srcs,ralb,rahb,declb,dechb,ra_centre,dec_centre
time.sleep(5)
srcs_bright = srcs[srcs[:,0]>=flux_thres]
allsrcs=srcs_bright
print "srcsbright:::::::",srcs_bright
time.sleep(5)
Ppath = '%s/pfiles_small'%Program_dir
###Ppath='/sharefs/hbkg/user/saina/pfiles_small/'
os.system("rm -r %s/pfiles_%s_load"%(Ppath,outpath[-14:-1]))
os.system("mkdir %s/pfiles_%s_load"%(Ppath,outpath[-14:-1]))###old ObsID = inpath[-17:-4]
###os.environ["PFILES"]="%s/pfiles_%s_load;/home/hxmt/hxmtsoft2/soft/install/x86_64-unknown-linux-gnu-libc2.12/syspfiles"%(Ppath,outpath[-14:-1])
os.environ["PFILES"]="%s/pfiles_%s_load;/hxmt/soft/Astro/heasoft-6.27.2/x86_64-pc-linux-gnu-libc2.17/syspfiles"%(Ppath,outpath[-14:-1])


#if not os.Wpath.exists('%s'%outfile5):
#	print 'reload cpha!'
#	import createpha_betax as createpha
os.chdir(outpath)
createpha.cpha(infile0,"%s0"%outfile5)
createpha.cpha(infile1,"%s1"%outfile5)
createpha.cpha(infile2,"%s2"%outfile5)
#else:
#	print 'cpha done!'

numb_srcs =  len(srcs_bright)
AllModels.clear()
temp= "powerlaw"
temp +="+hxmtpsf"
k=1
##############load data and fit###########################################
AllData("1:1 %s0.grp 2:2 %s1.grp 3:3 %s2.grp"%(outfile5,outfile5,outfile5))
print "alldata loaded 1" ,infilestr
print "Fit Begin"
time.sleep(5)
d1 = AllData(1)
d2 = AllData(2)
d3 = AllData(3)
AllModels+=temp
m1 = AllModels(1)
m2 = AllModels(2)
m3 = AllModels(3)
ra0=rau=[]
dec0=decu=[]
nor0=noru=[]
for src in srcs_bright:
	ra,dec,norm = src[1],src[2],src[0]
	m1.setPars({1:"0 -1",2:"0 0.1"})
	m1.setPars({3:" 1		  -1   "})
	m1.setPars({4:"60	-0.001"})
	m1.setPars({5:" %f		  -1   "%ra})
	m1.setPars({6:" %f		  -1   "%dec})
	m1.setPars({7:" %f		  1   "%int(norm/15)})
	m2.setPars({1:"0 -1",2:"0 0.1"})
	m2.setPars({4:"0 -0.001"})
	m3.setPars({1:"0 -1",2:"0 0.1"})
	m3.setPars({4:"-60  -0.001"})
	#AllModels.show()
	Fit.query='yes'
	Fit.perform()
	Fit.error("maximum 350 1.0 7")
	###Fit.show()
	if abs(m1(7).error[1]-m1(7).error[0])>0:
		snr_0=m1(7).values[0]*2/abs(m1(7).error[1]-m1(7).error[0])
		if (m1(7).values[0]>2)&(snr_0>2):
			ra0=np.append(ra0,m1(5).values[0])
			dec0=np.append(dec0,m1(6).values[0])
			nor0=np.append(nor0,m1(7).values[0])
		else:
			rau=np.append(rau,m1(5).values[0])
			decu=np.append(decu,m1(6).values[0])
			noru=np.append(noru,m1(7).values[0])
	else:
		rau=np.append(rau,m1(5).values[0])
		decu=np.append(decu,m1(6).values[0])
		noru=np.append(noru,m1(7).values[0])

######################################################################################
#AllModels.clear()
src_bright1=np.array([nor0,ra0,dec0])
src_bright0=src_bright1.T
numb_srcs0 =  len(src_bright0)
if len(ra0)>10:
	psrcs=[[]]*int(1+numb_srcs0/10)
	for i in range(len(psrcs)):
		psrcs[i]=src_bright0[i*10:(i+1)*10]
else:
	psrcs=[src_bright0]

ra1=[]
dec1=[]
nor1=[]
for j in psrcs:
	temp= "powerlaw"
	k = 0
	for i, src in enumerate(j):
		temp +="+hxmtpsf"
		k=k+1
		print "srcs_bright ObsID : ",k,temp
	
	number=k-1
	parnumb = (k-1)*5+7
	##############load data and fit###########################################
	print "Fit Begin"
	time.sleep(5)
	AllModels+=temp
	m1 = AllModels(1)
	m2 = AllModels(2)
	m3 = AllModels(3)
	m1.setPars({1:"0 -1",2:"0 0.1"})
	for i, src in enumerate(j):
		ra,dec,norm = src[1],src[2],src[0]
		m1.setPars({i*5+3:" 1		  -1   "})
		m1.setPars({i*5+4:"60	-0.001"})
		m1.setPars({i*5+5:" %f		  -1   "%ra})
		m1.setPars({i*5+6:" %f		  -1   "%dec})
		m1.setPars({i*5+7:" %f		  1   "%int(norm/15)})
		m2.setPars({1:"0 -1",2:"0 0.1"})
		m2.setPars({i*5+4:"0 -0.001"})
		m3.setPars({1:"0 -1",2:"0 0.1"})
		m3.setPars({i*5+4:"-60  -0.001"})
	
	#AllModels.show()
	Fit.query='yes'
	Fit.perform()
	Fit.error("maximum 350 1.0 1-%s"%(parnumb+1))
	###Fit.show()
	for i in range(7,parnumb+1,5):
		if abs(m1(i).error[1]-m1(i).error[0])>0:
			snr_1=m1(i).values[0]*2/abs(m1(i).error[1]-m1(i).error[0])
			if (m1(i).values[0]>2)&(snr_1>2):
				ra1=np.append(ra1,m1(i-2).values[0])
				dec1=np.append(dec1,m1(i-1).values[0])
				nor1=np.append(nor1,m1(i).values[0])
			else:
				rau=np.append(rau,m1(i-2).values[0])
				decu=np.append(decu,m1(i-1).values[0])
				noru=np.append(noru,m1(i).values[0])
		else:
			rau=np.append(rau,m1(i-2).values[0])
			decu=np.append(decu,m1(i-1).values[0])
			noru=np.append(noru,m1(i).values[0])
	
	ra0=ra1
	dec0=dec1
	nor0=nor1
	#AllModels.clear()

##################################################################################################
src_new=np.array([nor0,ra0,dec0])
srcs_bright=src_new.T
numb_srcs =  len(srcs_bright)
#AllModels.clear()
temp= "powerlaw"
k = 0
for i, src in enumerate(srcs_bright):
	temp +="+hxmtpsf"
	k=k+1
	print "srcs_bright ObsID : ",k,temp

number=k-1
add_ra=[]
add_dec=[]
addnum = 1
for i in range(addnum):
	temp +="+hxmtpsf"

temp+="+hxmtpsf"
parnumb = (k-1)*5+7
srclist = []
Lrarg = int(ralb-4)
Hrarg = int(rahb+4)
Ldecrg = int(declb-4)
Hdecrg = int(dechb+4)

atthd = pf.open(infile3)
atttb = atthd[1].data
atttm = atttb.field(0)

tbatt = atttb
ratp = tbatt.field(1)
dectp = tbatt.field(2)

#print "Fit Begin"
AllModels+=temp
m1 = AllModels(1)
m2 = AllModels(2)
m3 = AllModels(3)
m1.setPars({1:"0 -1",2:"0 0.1"})

for i, src in enumerate(srcs_bright):
	ra,dec,norm = src[1],src[2],src[0]
	m1.setPars({i*5+3:" 1		  -1   "})
	m1.setPars({i*5+4:"60	-0.001"})
	m1.setPars({i*5+5:" %f		  -1   "%ra})
	m1.setPars({i*5+6:" %f		  -1   "%dec})
	m1.setPars({i*5+7:" %f		  -1   "%int(norm/15)})
	m2.setPars({1:"0 -1",2:"0 0.1"})
	m2.setPars({i*5+4:"0 -0.001"})
	m3.setPars({1:"0 -1",2:"0 0.1"})
	m3.setPars({i*5+4:"-60  -0.001"})

for i in range(addnum):
	m1.setPars({i*5+5*k+3:" 1		  -1   "})
	m1.setPars({i*5+5*k+4:"60	-0.001"})
	m1.setPars({i*5+5*k+5:"%f,-1.0,%f,%f,%f,%f"%(ra_centre,Lrarg,Lrarg,Hrarg,Hrarg)})
	m1.setPars({i*5+5*k+6:"%f,-1.0,%f,%f,%f,%f"%(dec_centre,Ldecrg,Ldecrg,Hdecrg,Hdecrg)})
	m1.setPars({i*5+5*k+7:" %f		  -1   "%0.})
	m2.setPars({i*5+5*k+4:"0 -0.001"})
	m3.setPars({i*5+5*k+4:"-60  -0.001"})

m1.setPars({addnum*5+5*k+3:" 1		  -1   "})
m1.setPars({addnum*5+5*k+4:"60	-0.001"})
m1.setPars({addnum*5+5*k+5:"0	-1"})
m1.setPars({addnum*5+5*k+6:"0	-1"})
m1.setPars({addnum*5+5*k+7:"0	-1"})
m2.setPars({addnum*5+5*k+4:"0 -0.001"})
m3.setPars({addnum*5+5*k+4:"-60  -0.001"})

AllModels.show()
Fit.query='yes'
blist=ublist=[]
for i in range(7,parnumb+1,5):
	m1(i).frozen = False

Fit.perform()
Fit.error("maximum 350 1.0 1-%s"%(parnumb+1))
###Fit.show()

for i in range(7,parnumb+1,5):
	if abs(m1(i).error[1]-m1(i).error[0])>0:
		snr_2=m1(i).values[0]*2/abs(m1(i).error[1]-m1(i).error[0])
	else:
		snr_2=0
	if (m1(i).values[0]<=2)|(snr_2<=2):
		ublist=np.append(ublist,i)
		m1.setPars({i:"0 -1"})
	else:
		blist=np.append(blist,i)

if len(blist)>0:
	blist=blist.astype(int)

if len(ublist)>0:
	ublist=ublist.astype(int)

AllModels.show()
Fit.query='yes'
Fit.perform()

Plot.device = "%s/cps"%outfile2
Plot("data chi")
# Plot using Matplotlib:
fig = plt.figure(figsize = (20,10))
gs = gridspec.GridSpec(3, 1,bottom=0.05,top=0.96,left=0.05,right=0.95,hspace=0.05)
gs0 = gridspec.GridSpecFromSubplotSpec(2, 1,subplot_spec=gs[0], hspace=0.)
ax00 = plt.subplot(gs0[0])
ax01 = plt.subplot(gs0[1],sharex=ax00)
gs1 = gridspec.GridSpecFromSubplotSpec(2, 1,subplot_spec=gs[1], hspace=0.)
ax10 = plt.subplot(gs1[0],sharex=ax00)
ax11 = plt.subplot(gs1[1],sharex=ax00)
gs2 = gridspec.GridSpecFromSubplotSpec(2, 1,subplot_spec=gs[2], hspace=0.)
ax20 = plt.subplot(gs2[0],sharex=ax00)
ax21 = plt.subplot(gs2[1],sharex=ax00)
#axe[0].plot(d1.noticed, (np.array(d1.values)),'b-',d1.noticed, np.array(m1.folded(1)),'r-',lw=0.7)
#axe[0].errorbar(d1.noticed,np.array(d1.values),yerr=np.array(np.sqrt(d1.values)),fmt="o")
lc0=pf.open(infile0)
lc1=pf.open(infile1)
lc2=pf.open(infile2)
lcerror0 = lc0[1].data.field(2)
lcerror1 = lc1[1].data.field(2)
lcerror2 = lc2[1].data.field(2)
ax00.scatter(d1.noticed, (np.array(d1.values)),c='b',s=0.7)
ax00.plot(d1.noticed, np.array(m1.folded(1)),'r-',lw=0.7)
ax00.errorbar(d1.noticed, (np.array(d1.values)), yerr=np.array(lcerror0), zorder=0, fmt="none",marker="none",elinewidth=0.7)
ax01.plot(d1.noticed, (np.array(d1.values)-np.array(m1.folded(1))),lw=0.7)

ax10.scatter(d2.noticed, (np.array(d2.values)),c='b',s=0.7)
ax10.plot(d2.noticed, np.array(m2.folded(2)),'r-',lw=0.7)
ax10.errorbar(d2.noticed, (np.array(d2.values)), yerr=np.array(lcerror1), zorder=0, fmt="none",marker="none",elinewidth=0.7)
ax11.plot(d2.noticed, (np.array(d2.values)-np.array(m2.folded(2))),lw=0.7)

ax20.scatter(d3.noticed, (np.array(d3.values)),c='b',s=0.7)
ax20.plot(d3.noticed, np.array(m3.folded(3)),'r-',lw=0.7)
ax20.errorbar(d3.noticed, (np.array(d3.values)), yerr=np.array(lcerror2), zorder=0, fmt="none",marker="none",elinewidth=0.7)
ax21.plot(d3.noticed, (np.array(d3.values)-np.array(m3.folded(3))),lw=0.7)
ax21.set_xlabel('Times(s)',fontsize=18)
ax00.tick_params(labelbottom=False)
ax01.tick_params(labelbottom=False)
ax11.tick_params(labelbottom=False)
ax10.tick_params(labelbottom=False)
ax20.tick_params(labelbottom=False)
ax01.get_yticklabels()[-1].set_color('None')
ax11.get_yticklabels()[-1].set_color('None')
ax21.get_yticklabels()[-1].set_color('None')
plt.xlim(0,len(d1.noticed))
ax00.xaxis.grid(True)
ax01.xaxis.grid(True)
ax10.xaxis.grid(True)
ax11.xaxis.grid(True)
ax20.xaxis.grid(True)
ax21.xaxis.grid(True)
ax00.yaxis.grid(True)
ax01.yaxis.grid(True)
ax10.yaxis.grid(True)
ax11.yaxis.grid(True)
ax20.yaxis.grid(True)
ax21.yaxis.grid(True)
ax00.set_title('%s bright'%(outpath[-14:-1]))
axe = [ax00,ax10,ax20]
max1 = max(np.array(d1.values))
max2 = max(np.array(d2.values))
max3 = max(np.array(d3.values))
min1 = min(np.array(d1.values))
min2 = min(np.array(d2.values))
min3 = min(np.array(d3.values))
max1 = max1 + 18
max2 = max2 + 18
max3 = max3 + 18
maxt = [max1,max2,max3]
axe[0].set_ylim(min1-5,max1)
axe[1].set_ylim(min2-5,max2)
axe[2].set_ylim(min3-5,max3)
#plt.tight_layout(pad=0, h_pad=None, w_pad=None, rect=None)
plt.savefig(outpath+outpath[-14:-1]+'_bright.eps',dpi = 600)###png change to eps

#############################################################################################33
chi = Fit.statistic
dof = Fit.dof
rechi = chi/dof
for i in range(3,parnumb+1):
	m1(i).frozen = True

#############################################################################################  
print 'Begin to find new source....'
time.sleep(5)
#ratpp=ratp[::120]
#dectpp=dectp[::120]
file_newsource = outpath+outpath[-14:-1]+'_new.txt'
addtime=0
rechitemp=4
turn=0
for i in range(addnum):
	ratpp=ratp[::55];dectpp=dectp[::55]
	m1(i*5+5*k+5).frozen = False
	m1(i*5+5*k+6).frozen = False
	m1(i*5+5*k+7).frozen = False
	Fit.perform()
	print "****************** Start Addnum %s ! *********************"%(i+1)
	time.sleep(5)
	Fit.error("maximum 10 1.0 %s-%s"%((i*5+5*k+5),(5*k+7+5*i)))
	try:
		snr = 2*m1(i*5+5*k+7).values[0]/abs(m1(i*5+5*k+7).error[1]-m1(i*5+5*k+7).error[0])
	except:
		snr = 99999999
	print snr,m1(i*5+5*k+7).values[0]
	  
	x=0
	ratemp,dectemp,snrtemp,normtemp=0,0,0,0
	newra_try=[]
	newdec_try=[]
	radius=0.06
	for index,j in enumerate(ratpp):
		h = dectpp[index]
		if len(src_ra[((src_ra[:]-j)**2<=radius) & ((src_dec[:]-h)**2<=radius)])==0:  
			newra_try=np.r_[newra_try,j] 
			newdec_try=np.r_[newdec_try,h]
		print len(newra_try),len(newdec_try)  
	
	for index,j in enumerate(newra_try):
		h = newdec_try[index]	 
		m1.setPars({i*5+5*k+5:"%f,1.0,%f,%f,%f,%f"%(j,Lrarg,Lrarg,Hrarg,Hrarg)})
		m1.setPars({i*5+5*k+6:"%f,1.0,%f,%f,%f,%f"%(h,Ldecrg,Ldecrg,Hdecrg,Hdecrg)})
		Fit.perform()
		Fit.error("maximum 10 1.0 %s-%s"%(i*5+5*k+5,i*5+5*k+7))
		chi = Fit.statistic
		dof = Fit.dof
		rechi = chi/dof
	
		try:
			snr = 2*m1(i*5+5*k+7).values[0]/abs(m1(i*5+5*k+7).error[1]-m1(i*5+5*k+7).error[0])
			
		except:
			snr = 99999999
		 
		if ((rechi<rechitemp)&(snr!=99999999)&(m1(i*5+5*k+7).values[0]>0.9)):
			ratemp = m1(i*5+5*k+5).values[0]
			dectemp= m1(i*5+5*k+6).values[0]
			normtemp=m1(i*5+5*k+7).values[0]
			snrtemp= snr
			rechitemp=rechi
	#if ((snrtemp >=3.0)&(snrtemp < 99999)&(normtemp>0.6)):
	if normtemp >0:
		m1.setPars({i*5+5*k+5:"%s  1.0"%ratemp})
		m1.setPars({i*5+5*k+6:"%s  1.0"%dectemp})
		m1.setPars({i*5+5*k+7:"%s  1.0"%normtemp})
		Fit.perform()
		Fit.error("maximum 30 1.0 %s-%s"%(i*5+5*k+5,i*5+5*k+7))
		ra_chi=(m1(i*5+5*k+5).error[1]-m1(i*5+5*k+5).error[0])/2
		dec_chi=(m1(i*5+5*k+6).error[1]-m1(i*5+5*k+6).error[0])/2
		norm_err1=m1(i*5+5*k+7).error[0]
		norm_err2=m1(i*5+5*k+7).error[1]
		norm_err3=(m1(i*5+5*k+7).error[1]-m1(i*5+5*k+7).error[0])/2
		ra_values=m1(i*5+5*k+5).values[0]
		dec_values=m1(i*5+5*k+6).values[0]
		chi = Fit.statistic
		dof = Fit.dof
		rechi = chi/dof	
		with open(file_newsource,'a') as f:
				f.write('%s'%(i+1) +'\t''%s'%ra_values+'\t''%s'%ra_chi+'\t''%s'%dec_values+'\t''%s'%dec_chi+\
				'\t''%s'%m1(i*5+5*k+7).values[0]+'\t''%s'%norm_err1+'\t''%s'%norm_err2+'\t''%s'%norm_err3+\
				'\t''%s'%snrtemp+'\t''%s'%rechi+'\n')
		try:
			snr = 2*m1(i*5+5*k+7).values[0]/abs(m1(i*5+5*k+7).error[1]-m1(i*5+5*k+7).error[0])
		except:
			snr =99999999
		
		m1(i*5+5*k+5).frozen = True
		m1(i*5+5*k+6).frozen = True
		m1(i*5+5*k+7).frozen = True
		addtime = addtime+1
		if snr!=99999999:
			turn=turn+1
	else:
		#with open(file_newsource,'a') as f:
			#f.write('Not found New source!!!!!!!!!!!!!!!!!!'+'\n')
		print "##########################Do not found new source##############################"
		m1(i*5+5*k+5).frozen = True
		m1(i*5+5*k+6).frozen = True
		m1.setPars({i*5+5*k+7:"0  -1"})
		continue

###m1.show()
###################################################################################################
if rechi>50:
	print "reduced chi2 >20 the fit is bad."

for i in blist:
	m1(i).frozen = False

for i in range(addnum):
	m1(parnumb+(i+1)*5).frozen=False

Fit.perform()
Fit.error("maximum 350 1.0 1-%s"%(parnumb+1+addnum*5))
###Fit.show()
chi1 = Fit.statistic
dof1 = Fit.dof
notfitlist=[]
fitlist=[]
for i in blist:
	if (abs(m1(i).error[1]-m1(i).error[0])>0):
		snr_3=m1(i).values[0]*2/abs(m1(i).error[1]-m1(i).error[0])
	else:
		snr_3=0
	if (m1(i).values[0]<=2)|(snr_3<=2):
		notfitlist.append(i)
		rau=np.append(rau,m1(i-2).values[0])
		decu=np.append(decu,m1(i-1).values[0])
		noru=np.append(noru,m1(i).values[0])
		m1.setPars({i:"0 -1"})
	else:
		fitlist.append(i)

for i in range(3,parnumb+1+addnum*5):
	m1(i).frozen = True

if len(notfitlist)>0:
	for i in fitlist:
		m1(i).frozen = False
	
	for i in range(addnum):
		m1(parnumb+(i+1)*5).frozen=False
	
	Fit.perform()

chi0 = Fit.statistic
chi = Fit.statistic
dof = Fit.dof
'''for i in notfitlist:
	m1(i).frozen = False
	#m1(i-1).frozen = False
	#m1(i-2).frozen = False
	Fit.perform()
	Fit.error("maximum 50 1.0 %s"%i)
	tempchi = Fit.statistic
	ftest = Fit.ftest(tempchi,dof-3,chi,dof)
	if  ftest<0.1:
	#m1(i-1).frozen = True
	#m1(i-2).frozen = True
		m1(i).frozen = True
		Fit.perform()
		chi = Fit.statistic
	if  ftest>=0.1:
		m1(i).frozen = True
		continue'''

for i in range(1,parnumb+1):
	print>>parFile,m1(i).values[0],m1(i).error[0],m1(i).error[1]

dof1 = Fit.dof
print>>parFile,"chi2 and dof for the fit: ",chi,dof1
srclist = []
srcpot=[]
unsrclist=[]
for i in fitlist:
	srclist.append(i)

for j in range(parnumb+5,parnumb+1+addnum*5,5):
	unsrclist.append(j)

##############################################################################################
j=0
for i in srclist:
	j=j+1
	for m in range(1,parnumb+1+addnum*5):
		m1(m).frozen = True
	m1(2).frozen = False
	m1(i).frozen = False
	Fit.perform()
	Fit.error("maximum 50 1.0 %s-%s"%(i-2,i))
	fig,axe = plt.subplots(3,1,figsize=(16,9))
	tempmi = m1(i).values[0]
	fdtemp1 = np.array(m1.folded(1))
	fdtemp2 = np.array(m2.folded(2))
	fdtemp3 = np.array(m3.folded(3))
	axe[0].plot(d1.noticed, (np.array(d1.values)),'b-',d1.noticed, np.array(m1.folded(1)),'r-',lw=0.7)
	axe[1].plot(d2.noticed, (np.array(d2.values)),'b-',d2.noticed, np.array(m2.folded(2)),'r-',lw=0.7)
	axe[2].plot(d3.noticed, (np.array(d3.values)),'b-',d3.noticed, np.array(m3.folded(3)),'r-',lw=0.7)
	m1.setPars({i:"0 -1"})
	axe[0].plot(d1.noticed, np.array(m1.folded(1)),'g-',lw=0.7)
	axe[1].plot(d2.noticed, np.array(m2.folded(2)),'g-',lw=0.7)
	axe[2].plot(d3.noticed, np.array(m3.folded(3)),'g-',lw=0.7)
	axe[0].xaxis.grid(True)
	axe[1].xaxis.grid(True)
	axe[2].xaxis.grid(True)
	axe[0].set_title('%s  model%s  ra = %s dec= %s nor=%s'%(outpath[-14:-1],(i-2)/5,m1(i-2).values[0],m1(i-1).values[0],m1(i).values[0]))
	plt.savefig(outpath+outpath[-14:-1]+'_src%s.eps'%(j),dpi = 400)###png
	per_model=np.vstack((d1.noticed,np.array(d1.values),fdtemp1-np.array(m1.folded(1)),\
		d2.noticed,np.array(d2.values),fdtemp2-np.array(m2.folded(2)),d3.noticed,np.array(d3.values),\
		fdtemp3-np.array(m3.folded(3))))
	np.savetxt(outpath+"model_%s.dat"%(j),per_model.T,fmt="%s")
	m1.setPars({i:"%s -1"%tempmi})

for i in unsrclist:
	j=j+1
	for m in range(1,parnumb+1+addnum*5):
		m1(m).frozen = True
	m1(2).frozen = False
	m1(i).frozen = False
	Fit.perform()
	Fit.error("maximum 50 1.0 %s-%s"%(i-2,i))
	fig,axe = plt.subplots(3,1,figsize=(16,9))
	tempmi = m1(i).values[0]
	fdtemp1 = np.array(m1.folded(1))
	fdtemp2 = np.array(m2.folded(2))
	fdtemp3 = np.array(m3.folded(3))
	axe[0].plot(d1.noticed, (np.array(d1.values)),'b-',d1.noticed, np.array(m1.folded(1)),'r-',lw=0.7)
	axe[1].plot(d2.noticed, (np.array(d2.values)),'b-',d2.noticed, np.array(m2.folded(2)),'r-',lw=0.7)
	axe[2].plot(d3.noticed, (np.array(d3.values)),'b-',d3.noticed, np.array(m3.folded(3)),'r-',lw=0.7)
	m1.setPars({i:"0 -1"})
	axe[0].plot(d1.noticed, np.array(m1.folded(1)),'g-',lw=0.7)
	axe[1].plot(d2.noticed, np.array(m2.folded(2)),'g-',lw=0.7)
	axe[2].plot(d3.noticed, np.array(m3.folded(3)),'g-',lw=0.7)
	axe[0].xaxis.grid(True)
	axe[1].xaxis.grid(True)
	axe[2].xaxis.grid(True)
	axe[0].set_title('%s  model%s  ra = %s dec= %s nor=%s'%(outpath[-14:-1],(i-2)/5,m1(i-2).values[0],m1(i-1).values[0],m1(i).values[0]))
	plt.savefig(outpath+outpath[-14:-1]+'_unsrc%s.eps'%(j),dpi = 600)###eps
	per_model=np.vstack((d1.noticed,np.array(d1.values),fdtemp1-np.array(m1.folded(1)),\
		d2.noticed,np.array(d2.values),fdtemp2-np.array(m2.folded(2)),d3.noticed,np.array(d3.values),\
		fdtemp3-np.array(m3.folded(3))))
	np.savetxt(outpath+"ns_model_%s.dat"%(j),per_model.T,fmt="%s")
	m1.setPars({i:"%s -1"%tempmi})

###############################output to file and show picture####################################
Plot.device = "%s/cps"%outfile2
Plot("data chi")
# Plot using Matplotlib:
fig = plt.figure(figsize = (20,10))
gs = gridspec.GridSpec(3, 1,bottom=0.05,top=0.96,left=0.05,right=0.95,hspace=0.05)
gs0 = gridspec.GridSpecFromSubplotSpec(2, 1,subplot_spec=gs[0], hspace=0.)
ax00 = plt.subplot(gs0[0])
ax01 = plt.subplot(gs0[1],sharex=ax00)
gs1 = gridspec.GridSpecFromSubplotSpec(2, 1,subplot_spec=gs[1], hspace=0.)
ax10 = plt.subplot(gs1[0],sharex=ax00)
ax11 = plt.subplot(gs1[1],sharex=ax00)
gs2 = gridspec.GridSpecFromSubplotSpec(2, 1,subplot_spec=gs[2], hspace=0.)
ax20 = plt.subplot(gs2[0],sharex=ax00)
ax21 = plt.subplot(gs2[1],sharex=ax00)
#axe[0].plot(d1.noticed, (np.array(d1.values)),'b-',d1.noticed, np.array(m1.folded(1)),'r-',lw=0.7)
#axe[0].errorbar(d1.noticed,np.array(d1.values),yerr=np.array(np.sqrt(d1.values)),fmt="o")
ax00.scatter(d1.noticed, (np.array(d1.values)),c='b',s=0.7)
ax00.plot(d1.noticed, np.array(m1.folded(1)),'r-',lw=0.7)
ax00.errorbar(d1.noticed, (np.array(d1.values)), yerr=np.array(lcerror0), zorder=0, fmt="none",marker="none",elinewidth=0.7)
ax01.plot(d1.noticed, (np.array(d1.values)-np.array(m1.folded(1))),lw=0.7)

ax10.scatter(d2.noticed, (np.array(d2.values)),c='b',s=0.7)
ax10.plot(d2.noticed, np.array(m2.folded(2)),'r-',lw=0.7)
ax10.errorbar(d2.noticed, (np.array(d2.values)), yerr=np.array(lcerror1), zorder=0, fmt="none",marker="none",elinewidth=0.7)
ax11.plot(d2.noticed, (np.array(d2.values)-np.array(m2.folded(2))),lw=0.7)

ax20.scatter(d3.noticed, (np.array(d3.values)),c='b',s=0.7)
ax20.plot(d3.noticed, np.array(m3.folded(3)),'r-',lw=0.7)
ax20.errorbar(d3.noticed, (np.array(d3.values)), yerr=np.array(lcerror2), zorder=0, fmt="none",marker="none",elinewidth=0.7)
ax21.plot(d3.noticed, (np.array(d3.values)-np.array(m3.folded(3))),lw=0.7)
ax21.set_xlabel('Times(s)',fontsize=18)
ax00.tick_params(labelbottom=False)
ax01.tick_params(labelbottom=False)
ax11.tick_params(labelbottom=False)
ax10.tick_params(labelbottom=False)
ax20.tick_params(labelbottom=False)
ax01.get_yticklabels()[-1].set_color('None')
ax11.get_yticklabels()[-1].set_color('None')
ax21.get_yticklabels()[-1].set_color('None')
plt.xlim(0,len(d1.noticed))
ax00.xaxis.grid(True)
ax01.xaxis.grid(True)
ax10.xaxis.grid(True)
ax11.xaxis.grid(True)
ax20.xaxis.grid(True)
ax21.xaxis.grid(True)
ax00.yaxis.grid(True)
ax01.yaxis.grid(True)
ax10.yaxis.grid(True)
ax11.yaxis.grid(True)
ax20.yaxis.grid(True)
ax21.yaxis.grid(True)
axe = [ax00,ax10,ax20]
max1 = max(np.array(d1.values))
max2 = max(np.array(d2.values))
max3 = max(np.array(d3.values))
min1 = min(np.array(d1.values))
min2 = min(np.array(d2.values))
min3 = min(np.array(d3.values))
max1 = max1 + 18
max2 = max2 + 18
max3 = max3 + 18
maxt = [max1,max2,max3]
axe[0].set_ylim(min1-5,max1)
axe[1].set_ylim(min2-5,max2)
axe[2].set_ylim(min3-5,max3)
ax00.set_title('%s bright&new'%(outpath[-14:-1]))
plt.savefig(outpath+outpath[-14:-1]+'_bn.eps',dpi = 600)###png2eps
print "step1 end"
#######################################################################################################
allsrc=np.r_[srclist,unsrclist].astype(int)
for i in allsrc:
	m1(i).frozen=False

Fit.perform()
#re_chi = chi/dof1
print chi/dof1,rechi
max_chi = 3
if rechi>max_chi:
	max_chi = rechi+0.01

#chi_2 = np.vstack((chi,dof1,rechi))
#np.savetxt("./result/chi2.dat",chi_2.T,fmt="%s")
for i in allsrc:
	Fit.error("max %s 1.0 %s"%(max_chi,i))

data_all = np.vstack((d1.noticed,np.array(d1.values),np.array(m1.folded(1)),\
np.array(d1.values)-np.array(m1.folded(1)),d2.noticed,np.array(d2.values),np.array(m2.folded(2)),\
np.array(d2.values)-np.array(m2.folded(2)),d3.noticed,np.array(d3.values),np.array(m3.folded(3)),\
np.array(d3.values)-np.array(m3.folded(3))))
np.savetxt(outpath+outpath[-14:-1]+"_data_3box.dat",data_all.T,fmt="%s")
chi_2 = np.vstack((chi,dof1,rechi,m1(2).values[0],m1(2).error[0],m1(2).error[1],m2(2).values[0],\
m2(2).error[0],m2(2).error[1],m3(2).values[0],m3(2).error[0],m3(2).error[1]))
np.savetxt(outpath+outpath[-14:-1]+"_chi2.dat",chi_2.T,fmt="%s")
print "step2 end"
############################################################################################################
values=[]
ra=[]
dec=[]
error1=[]
error2=[]
error3=[]
snr=[]
names=[]
number1=[]
number2=[]
if parnumb<3:
	print "No bright source in src_map"
else:
	for i in srclist:
		if m1(i).values[0]<2:
			rau=np.append(rau,m1(i-2).values[0])
			decu=np.append(decu,m1(i-1).values[0])
			noru=np.append(noru,m1(i).values[0])
			continue
		if (m1(i).error[1])-(m1(i).error[0])>0:
			snr1=m1(i).values[0]*2/abs(m1(i).error[1]-m1(i).error[0])
		else:
			rau=np.append(rau,m1(i-2).values[0])
			decu=np.append(decu,m1(i-1).values[0])
			noru=np.append(noru,m1(i).values[0])
			continue
		values=np.append(values,m1(i).values[0])
		ra=np.append(ra,m1(i-2).values[0])
		dec=np.append(dec,m1(i-1).values[0])
		error1=np.append(error1,m1(i).error[0])
		error2=np.append(error2,m1(i).error[1])
		error3=np.append(error3,((m1(i).error[1]-m1(i).error[0])/2))
		snr=np.append(snr,snr1)

if ra==[]:
	print "No bright source in hxmt_map"
else:
	f1=open(outpath+outpath[-14:-1]+"_source_bright.txt","w+")
	f2=open(outpath+outpath[-14:-1]+"_source_moreinfo.txt","w+")
	sn=len(ra)
	for i in range(sn):
		for j in range(len(src_ra)):
			if (src_ra[j]>=ra[i]-0.000001)&(src_ra[j]<=ra[i]+0.000001)&\
			(src_dec[j]>=dec[i]-0.000001)&(src_dec[j]<=dec[i]+0.000001):
				print j+1,src_name[j],src_ra[j],src_dec[j]
				names=np.append(names,src_name[j])
				number1=np.append(number1,j+1)
				number2=np.append(number2,i+1)
				newstr2=[str(int(i+1)),str(int(j+1)),src_name[j],str(ra[i]),str(dec[i]),\
				str(values[i]),str(error1[i]),str(error2[i]),str(error3[i]),str(snr[i])]
				newstr1=[str(int(i+1)),str(int(j+1)),src_name[j],str(ra[i]),str(dec[i]),\
				str(values[i]),str(error3[i]),str(snr[i])]
				filestr1='\t'.join(newstr1)
				filestr2='\t'.join(newstr2)
				f2.writelines(filestr2+'\n')
				f1.writelines(filestr1+'\n')
	
	f1.close()
	f2.close()

############################################################################################
res0 = np.array(d1.values)-np.array(m1.folded(1))
res1 = np.array(d2.values)-np.array(m2.folded(2))
res2 = np.array(d3.values)-np.array(m3.folded(3))
lc_all0 = lc0
lc_all1 = lc1
lc_all2 = lc2
lc_data0 = lc_all0[1].data
lc_gti = lc_all0[2]
lc_time0 = lc_data0['Time']
lc_ct0 = lc_data0['Counts']
lc_ct_e0 = lc_data0['Stat_err']
lc_data1 = lc_all1[1].data
lc_time1 = lc_data1.field(0)
lc_ct1 = lc_data1.field(1)
lc_ct_e1 = lc_data1.field(2)
lc_data2 = lc_all2[1].data
lc_time2 = lc_data2.field(0)
lc_ct2 = lc_data2.field(1)
lc_ct_e2 = lc_data2.field(2)
col10 = pf.Column(name = 'Time', format='D',array = lc_time0)
col20 = pf.Column(name = 'Counts', format='D',array = res0)
col30 = pf.Column(name = 'Stat_err', format='D',array = lc_ct_e0)
cols0 = pf.ColDefs([col10,col20,col30])
tbhdu0 = pf.BinTableHDU.from_columns(cols0)
prihdu = pf.PrimaryHDU()
gtihdu = lc_gti
tbhdulist0 = pf.HDUList([prihdu,tbhdu0,gtihdu])
tbhdulist0.writeto(outpath+"me_lc_box0_small_res.fits",clobber=True)

col11 = pf.Column(name = 'Time', format='D',array = lc_time1)
col21 = pf.Column(name = 'Counts', format='D',array = res1)
col31 = pf.Column(name = 'Stat_err', format='D',array = lc_ct_e1)
cols1 = pf.ColDefs([col11,col21,col31])
tbhdu1 = pf.BinTableHDU.from_columns(cols1)
tbhdulist1 = pf.HDUList([prihdu,tbhdu1,gtihdu])
tbhdulist1.writeto(outpath+"me_lc_box1_small_res.fits",clobber=True)

col12 = pf.Column(name = 'Time', format='D',array = lc_time2)
col22 = pf.Column(name = 'Counts', format='D',array = res2)
col32 = pf.Column(name = 'Stat_err', format='D',array = lc_ct_e2)
cols2 = pf.ColDefs([col12,col22,col32])
tbhdu2 = pf.BinTableHDU.from_columns(cols2)
tbhdulist2 = pf.HDUList([prihdu,tbhdu2,gtihdu])
tbhdulist2.writeto(outpath+"me_lc_box2_small_res.fits",clobber=True)
print "res: finish"
##########################################################################################
'''resfile0=outpath+'me_lc_box0_small_res.fits'
resfile1=outpath+'me_lc_box1_small_res.fits'
resfile2=outpath+'me_lc_box2_small_res.fits'
createpha.cpha(resfile0,"%s0"%outfile5)
createpha.cpha(resfile1,"%s1"%outfile5)
createpha.cpha(resfile2,"%s2"%outfile5)

AllModels.clear()
temp= "powerlaw"
temp +="+hxmtpsf"
k=1
AllData("1:1 %s0.grp 2:2 %s1.grp 3:3 %s2.grp"%(outfile5,outfile5,outfile5))
print "Fit UB Begin"
time.sleep(2)
d1 = AllData(1)
d2 = AllData(2)
d3 = AllData(3)
AllModels+=temp
m1 = AllModels(1)
m2 = AllModels(2)
m3 = AllModels(3)'''
for i in range(1,parnumb+addnum*5+5+1):
	m1(i).frozen=True

aulist0=np.array([noru,rau,decu])
aulist1=aulist0.T[np.lexsort(aulist0[1:2,:])]
aulist1[:,0]=0
aulist=np.unique(aulist1,axis=0)
numb_srcs=len(aulist)
f4=open(outpath+outpath[-14:-1]+"_source_ub.txt","w+")
f3=open(outpath+outpath[-14:-1]+"_ub_moreinfo.txt","w+")
values1=ra1=dec1=error11=error21=error31=snr1=names1=number1=number21=[]
chi_2=[[]]*numb_srcs
m=0
eflist=[]
#arslist=[]
bkg1=m1(2).values[0];bkg2=m2(2).values[0];bkg3=m3(2).values[0]
m1.setPars({2:"%f,1.0,%f,%f,%f,%f"%(bkg1,bkg1-0.001,bkg1-0.001,bkg1+0.01,bkg1+0.01)})
m2.setPars({2:"%f,1.0,%f,%f,%f,%f"%(bkg2,bkg2-0.001,bkg2-0.001,bkg2+0.01,bkg2+0.01)})
m3.setPars({2:"%f,1.0,%f,%f,%f,%f"%(bkg3,bkg3-0.001,bkg3-0.001,bkg3+0.01,bkg3+0.01)})

for i,src in enumerate(aulist):
	ra,dec,norm = src[1],src[2],src[0]
	m1.setPars({addnum*5+5*k+5:" %f		  -1   "%ra})
	m1.setPars({addnum*5+5*k+6:" %f		  -1   "%dec})
	m1.setPars({addnum*5+5*k+7:"%f,1,-10000,-10000,10000,10000"%0.})
	#m1.show()
	l=addnum*5+5*k+7
	#m1(l).frozen = False
	Fit.query='yes'
	Fit.perform()
	chi = Fit.statistic
	dof = Fit.dof
	rechi = chi/dof
	if rechi>50:
		print "reduced chi2 >20 the fit is bad."
	
	Fit.error("maximum 50 1.0 %s"%l)
	###Fit.show()
	chi1 = Fit.statistic
	dof1 = Fit.dof
	if abs(m1(l).error[1])-(m1(l).error[0])>0:
		snru=abs(m1(l).values[0]*2/abs(m1(l).error[1]-m1(l).error[0]))
	else:
		continue
	
	efar=sum(np.array(m1.folded(1))-m1(2).values[0]+np.array(m2.folded(2))-m2(2).values[0]+np.array(m3.folded(3))-m3(2).values[0])/(m1(l).values[0])
	eflist=np.append(eflist,efar)
	print '-------------------------------',i,efar,'-------------------------------'
	#if efar<=0.5:
		#arslist=np.append(arlist,i)
		#continue
	chi_2[m] = chi,dof,rechi
	m=m+1
	values1=m1(l).values[0]
	ra1=m1(l-2).values[0]
	dec1=m1(l-1).values[0]
	error1=m1(l).error[0]
	error2=m1(l).error[1]
	error3=((m1(l).error[1]-m1(l).error[0])/2)
	snr1=snru
	data_all1 = np.vstack((d1.noticed,np.array(d1.values),np.array(m1.folded(1)),np.array(d1.values)-np.array(m1.folded(1)),\
		d2.noticed,np.array(d2.values),np.array(m2.folded(2)),np.array(d2.values)-np.array(m2.folded(2)),\
		d3.noticed,np.array(d3.values),np.array(m3.folded(3)),np.array(d3.values)-np.array(m3.folded(3))))
	np.savetxt(outpath+"ub_model_%s.dat"%(m),data_all1.T,fmt="%s")
	#print>>parFile,"chi2 and dof for the fit: ",chi,dof
	#Xset.closeLog()
	#fig,axe = plt.subplots(3,1,figsize=(16,9))
	#axe[0].plot(d1.noticed, (np.array(d1.values)),'b-',d1.noticed, np.array(m1.folded(1)),'r-',lw=0.7)
	#axe[1].plot(d2.noticed, (np.array(d2.values)),'b-',d2.noticed, np.array(m2.folded(2)),'r-',lw=0.7)
	#axe[2].plot(d3.noticed, (np.array(d3.values)),'b-',d3.noticed, np.array(m3.folded(3)),'r-',lw=0.7)
	#axe[0].xaxis.grid(True)
	#axe[1].xaxis.grid(True)
	#axe[2].xaxis.grid(True)
	#plt.show()
	for j in range(len(src_ra)):
		if (src_ra[j]>=ra1-0.000001)&(src_ra[j]<=ra1+0.000001)&(src_dec[j]>=dec1-0.000001)&(src_dec[j]<=dec1+0.000001):
			print j+1,src_name[j],src_ra[j],src_dec[j]
			names1=np.append(names1,src_name[j])
			number1=np.append(number1,j+1)
			number21=np.append(number21,m)
			newstr1=[str(int(m)),str(int(j+1)),src_name[j],str(ra1),str(dec1),str(m1(l).values[0]),\
				str(m1(l).error[0]),str(m1(l).error[1]),str((m1(l).error[1]-m1(l).error[0])/2),str(snr1)]
			filestr1='\t'.join(newstr1)
			f3.writelines(filestr1+'\n')
			newstr2=[str(int(m)),str(int(j+1)),src_name[j],str(ra1),str(dec1),str(m1(l).values[0]),\
				str((m1(l).error[1]-m1(l).error[0])/2),str(snr1)]
			filestr2='\t'.join(newstr2)
			f4.writelines(filestr2+'\n')

#print aulist[arslist]
f3.close()
f4.close()
np.savetxt(outpath+outpath[-14:-1]+"_eflist.txt",eflist,fmt="%s")
np.savetxt(outpath+outpath[-14:-1]+"_ub_chi2.dat",chi_2,fmt="%s")
parFile.close()
srcFile.close()
plt.close('all')
os.system("rm -r %s/pfiles_%s_load"%(Ppath,outpath[-14:-1]))
sys.path.append('%s/lcfit'%Program_dir)
###sys.path.append('/sharefs/hbkg/user/saina/tar/P0201013/')
import src_web_auto
src_web_auto.dataout(outpath[:-20],outpath[-14:-1])
time2 = time.time()
program_tree = '/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN'
with open(program_tree + '/gps_mission/lcfit_success.txt', 'a+') as f:
	print >> f, outpath[-14:-1]
tmptm = str(time.strftime('%Y.%m.%d-%H:%M',time.localtime(time.time())))
with open(program_tree + '/gps_mission/lcfit_success.txt', 'a+') as f:
	print >> f, (outpath[-14:-1]+ '\tlc fit finished!  %s\n' % tmptm)
print 'fit finished ',time2
print "time expended: ",(time2-time1)/3600.0,' hours'

import project_rdata
import src_map_auto
project_rdata.dataout(outpath[:-20], outpath[-14:-1])
src_map_auto.src_map(outpath[:-20], outpath[-14:-1])
