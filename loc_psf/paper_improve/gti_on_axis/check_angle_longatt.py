import matplotlib
matplotlib.use("agg")
from xspec import *
from math import *
import numpy as np
from astropy.io import fits as pf
from readxml import *
import time
import sys
import Quat as quat
#from testroll import *
import testroll
from xml.etree.ElementTree import ElementTree,Element
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

'''
To calibrate PSF with Crab in FOV alpha<0.8, beta<0.8.
No need to distinguish 3 box.
#  Crab, ra = 83.633, dec = 22.014
    H 1730-333, ra = 263.353, dec = -33.3888
    GX 354-0, ra = 262.991, dec = -33.834
1. calculate crt att with lc time
2. use Crab ra/dec to get alpha/beta list
3. return lc time dex
4. merge time for use 
'''


if len(sys.argv)<2:
    print("Need the config file!")
    cfg = "config.xml"
    #sys.exit(1)
else:
    cfg = sys.argv[1]

year = sys.argv[2]

if len(sys.argv)==4:
    select_box_dex = int(sys.argv[3])
else:
    select_box_dex = 0

#---------------make config------------------------#
def gen_config(att_list,lc_list,fits_dir,xmlpath='./config_me.xml',instru='LE'):
    exec_name = './command_'+instru.lower()+'_'+time.strftime("%y%m%d")+'.sh'

    tree = ElementTree()
    tree.parse(xmlpath)
    pathnodes = tree.findall("PATH/outpath")
    outnodes = tree.findall("Outfile_Name/outfilename")
    filenodes = tree.findall("Infile_Name/infilename/infilenamelist")
    instrunodes = tree.findall("Inst_Info/Instrument")
    outfile = fits_dir + "/config_%s.xml"%(instru.lower())
    print(outfile)
    filenodes[0].text = "\n  " + lc_list[0]+" \n"
    filenodes[1].text = "\n  " + lc_list[1]+" \n"
    filenodes[2].text = "\n  " + lc_list[2]+" \n"
    filenodes[3].text = '\n  ' + att_list[0] + ' \n'
    pathnodes[0].text = '\n  ' + fits_dir + '  \n'
    instrunodes[0].text = '\n ' +instru + ' \n'
    outnodes[0].text = '\n '+'/GAL_he_small \n'

    tree.write(outfile, encoding="utf-8", xml_declaration=True)
    with open(exec_name,'a+')as f:
        try:
            print('sh ./sub_task_psf.sh %s'%(outfile),file=f)
        except:
            print('cannot write down command!!!')

#---------------------load file-  -------------#
readcfg = loadDom(cfg)
infilestr = readcfg.getTagText("infilenamelist")
inpathstr = readcfg.getTagText("inpath")
outpathstr = readcfg.getTagText("outpath")
instr = readcfg.getTagText("Instrument")
#print inpathstr,infilestr
inpathstr = inpathstr.strip()
outpathstr = outpathstr.strip()
evtfilestr0 = infilestr.split()[0]
evtfilestr1 = infilestr.split()[1]
evtfilestr2 = infilestr.split()[2]
infilestr2 = infilestr.split()[-1]
infile = (inpathstr+infilestr2)
evtfile0 = (inpathstr+evtfilestr0)
evtfile1 = (inpathstr+evtfilestr1)
evtfile2 = (inpathstr+evtfilestr2)
evtfile_list = [evtfile0,evtfile1,evtfile2]
dir_front = outpathstr

evtfile = evtfile_list[select_box_dex]

print ("the scanning pointing file : ",infile)
print("the lc file :", evtfile_list)
instr = instr.strip()
instr = instr.split()[0]
instr = instr#.encode()

if instr == "HE":
    Talfa_bound=5.7 #means long_bound
    Tbeta_bound=1.1 #means short_bound
    dcm_b_f = np.array([[0,1,0],[0,0,1],[1,0,0]])
elif instr == "ME":
    Talfa_bound=4.0 #means long_bound
    Tbeta_bound=1.0 #means short_bound
    dcm_b_f = np.array([[0,0,-1],[0,1,0],[1,0,0]])
elif instr == "LE":
    Talfa_bound=6.0 #means long_bound
    Tbeta_bound=1.6 #means short_bound
    dcm_b_f = np.array([[0,0,1],[0,1,0],[-1,0,0]])
else:
    print("The instrment in config file is wrong.")
    sys.exit(1)

#---------------------Crab ra/dec--------------#
ra =  83.633
dec = 22.014

alpha_lim = 0.8
beta_lim = 0.8

box_roll = [60,0,-60]
roll=box_roll[select_box_dex]/180.0*np.pi
roll_indx = int(-box_roll[select_box_dex]/60+1)

#--------------------read fits-----------------#
atthd = pf.open(infile)
lchd = pf.open(evtfile)
atdt = atthd[3].data
lcdt = lchd[1].data
qtime = atdt.field(0)
lctime = lcdt.field(0)
lccounts = lcdt.field(1)
lcerr = lcdt.field(2)
tstart = lctime[0]
tstop = lctime[-1]

#----------------quat with lc time------------#
q1_list = atdt.field(1)[np.in1d(qtime,lctime)]
q2_list = atdt.field(2)[np.in1d(qtime,lctime)]
q3_list = atdt.field(3)[np.in1d(qtime,lctime)]


#-------------calculate alpha/beta------------#
mtx_list = []
alpha_list = []
beta_list = []
for i in range(0,q1_list.shape[0]):
    quat1 = [q1_list[i],q2_list[i],q3_list[i]]
    mtx_list.append(testroll.quat2mtx(quat1))

psai = 0
theta = 0
phi = 0
rot_quat =  quat.Quat((psai,theta,phi))
x,y,z = testroll.sph2cart(ra*np.pi/180,dec*np.pi/180,1)
xyz = [x,y,z]
dcm_b_f_rot = rot_quat.transform.T


# print('length: lighr curvefile, Quat file >>>',len(lctime),len(q1_list))
#
# atthd.close()
# lchd.close()
# del atthd,lchd


mtx_list= []
accept_dex = []
accept_time = []
alpha_b0 = []
beta_b0 = []
alpha_b1 = []
beta_b1 = []
alpha_b2 = []
beta_b2 = []

time_len = 100

delta_alfa_box = np.zeros(3)
delta_beta_box = np.zeros(3)

for i in range(0,time_len):
    quat1 = [q1_list[i],q2_list[i],q3_list[i]]
    mtx = testroll.quat2mtx(quat1)
    mtx_list.append(mtx)
    delta_alfa0, delta_beta0, xr, yr, zr = testroll.quat2delta2_rot(mtx, dcm_b_f, xyz, dcm_b_f_rot)
    ###cosz = np.sqrt(zr ** 2 / (xr ** 2 + yr ** 2 + zr ** 2))
    #if delta_alfa>alpha_lim and delta_beta>beta_lim and delta_alfa<1.0 and delta_beta<1.0:
    for select_box_dex in range(3):
        roll = box_roll[select_box_dex] / 180.0 * np.pi
        delta_alfa_box[select_box_dex] = delta_alfa0 * cos(roll) - delta_beta0 * sin(roll)
        delta_beta_box[select_box_dex] = delta_alfa0 * sin(roll) + delta_beta0 * cos(roll)
    accept_dex.append(i)
    accept_time.append(lctime[i])
    alpha_b0.append(delta_alfa_box[0])
    beta_b0.append(delta_beta_box[0])
    alpha_b1.append(delta_alfa_box[1])
    beta_b1.append(delta_beta_box[1])
    alpha_b2.append(delta_alfa_box[2])
    beta_b2.append(delta_beta_box[2])

tot_list = np.zeros((10,time_len))
tot_list[0,:] = np.array(alpha_b0)
tot_list[1,:] = np.array(beta_b0)
tot_list[2,:] = np.array(alpha_b1)
tot_list[3,:] = np.array(beta_b1)
tot_list[4,:] = np.array(alpha_b2)
tot_list[5,:] = np.array(beta_b2)
tot_list[6,:] = np.array(q1_list[0:time_len])
tot_list[7,:] = np.array(q2_list[0:time_len])
tot_list[8,:] = np.array(q3_list[0:time_len])
tot_list[9,:] = np.array(accept_time)
print(accept_time)
print(accept_time[65])
print(accept_time[80])


def stand_psfmodel(delta_alfa0, delta_beta0):
    PI = 3.14159265358979323846
    alfa_bound = Talfa_bound  ###np.abs(Talfa_bound + 0)#abound[roll_indx]
    beta_bound = Tbeta_bound  ###np.abs(Tbeta_bound + 0)#bbound[roll_indx]
    delta_alfa = np.abs(delta_alfa0)
    delta_beta = np.abs(delta_beta0)
    condi = np.logical_and(delta_alfa < alfa_bound, delta_beta < beta_bound)
    print(condi.shape)
    factor = np.where(condi, (1.0 - np.tan(delta_alfa / 180.0 * PI) / np.tan(alfa_bound / 180.0 * PI)) * \
                      (1.0 - np.tan(delta_beta / 180.0 * PI) / np.tan(beta_bound / 180.0 * PI)),
                      np.zeros_like(delta_beta0))
    factor = factor / (np.sqrt(
        np.tan(delta_alfa / 180.0 * PI) * np.tan(delta_alfa / 180.0 * PI) + np.tan(delta_beta / 180.0 * PI) * np.tan(
            delta_beta / 180 * PI) + 1))
    return factor

tot_list = np.transpose(tot_list)
np.savetxt('angle_100sAtt.txt',tot_list)
tot_list = np.transpose(tot_list)


norm = lccounts[0:time_len]
yerr = lcerr[0:time_len]
color_tem = 'black'
fig = plt.figure(figsize=(18, 10))  # 10,20
gs = gridspec.GridSpec(6, 3)
ax1 = plt.subplot(gs[0:4, :])
###ax2 = plt.subplot(gs[3:5, :], sharex=ax1)
ax3 = plt.subplot(gs[4:, :], sharex=ax1)
gs.update(hspace=0)
x_values = np.arange(0, time_len)
ax1.errorbar(x_values, norm, yerr, label='lc',fmt='_',ecolor=color_tem,mfc=color_tem, mfcalt=color_tem,mec=color_tem
             ,markersize=6,elinewidth=1,capsize=0)


ave_cts =  np.sum(norm)/len(norm)

for i in range(3):

    alpha = tot_list[int(i*2)+0,:]
    beta = tot_list[int(i*2)+1,:]
    psf = stand_psfmodel(alpha,beta)*ave_cts
    print(len(alpha),len(beta),len(norm),len(psf),len(yerr))
    if i==0:
        sigma = (norm-psf)/yerr
    ax1.plot(x_values, psf, '.', label='box%s'%str(i))

ax3.errorbar(x_values, sigma, 1, label='sigma', fmt='_', ecolor=color_tem, mfc=color_tem, mfcalt=color_tem,
             mec=color_tem
             , markersize=6, elinewidth=1, capsize=0)
ax1.legend()
###ax2.legend()
ax3.legend()

'''
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 1, 1)  # , projection='3d')
plt.plot(x_values, alpha,label='alpha')
plt.plot(x_values, beta,label='beta')
plt.errorbar(x_values, norm, yerr,label='lc')
plt.plot(x_values,psf,label='psf')
plt.legend()
'''

fig.savefig('angle_longatt.png')

# fig = plt.figure(figsize=plt.figaspect(0.5))
# ax = fig.add_subplot(1, 1, 1)  # , projection='3d')
# plt.plot(np.array(accept_alpha), np.array(accept_beta),"*")
# fig.savefig('%s_check_angle.png'%(instr))