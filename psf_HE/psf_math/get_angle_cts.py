import matplotlib
matplotlib.use("agg")
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

# if len(sys.argv)==4:
#     select_box_dex = int(sys.argv[3])
# else:
#     select_box_dex = 1

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

def dcm_roll_theta(theta):
    theta = np.deg2rad(theta)
    dcm = np.zeros((3,3))#np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
    dcm[0,0] = np.cos(theta)
    dcm[0,1] = -np.sin(theta)
    dcm[1,0] = np.sin(theta)
    dcm[1,1] = np.cos(theta)
    dcm[2,2] = 1
    return dcm

#---------------------Crab ra/dec--------------#
ra =  83.633
dec = 22.014

# alpha_uplim = [-0.1,0.8,0.8]
# alpha_downlim = [-0.8,0.1,0.1]
# beta_uplim = [-0.1,-0.1,-0.1]
# beta_downlim = [-0.8,-0.8,-0.8]
alpha_uplim = [Talfa_bound,Talfa_bound,Talfa_bound]
alpha_downlim = [-Talfa_bound,-Talfa_bound,-Talfa_bound]
beta_uplim = [Tbeta_bound,Tbeta_bound,Tbeta_bound]
beta_downlim = [-Tbeta_bound,-Tbeta_bound,-Tbeta_bound]

box_roll = [60,0,-60]

###roll_indx = int(-box_roll[select_box_dex]/60+1)



lc_list = []


for select_box_dex in range(3):
    evtfile = evtfile_list[select_box_dex]
    # --------------------read fits-----------------#
    atthd = pf.open(infile)
    lchd = pf.open(evtfile)
    atdt = atthd[3].data
    lcdt = lchd[1].data
    qtime = atdt.field(0)
    lctime = lcdt.field(0)
    tstart = lctime[0]
    tstop = lctime[-1]

    # ----------------quat with lc time------------#
    q1_list = atdt.field(1)[np.in1d(qtime, lctime)]
    q2_list = atdt.field(2)[np.in1d(qtime, lctime)]
    q3_list = atdt.field(3)[np.in1d(qtime, lctime)]

    # -------------calculate alpha/beta------------#
    mtx_list = []
    alpha_list = []
    beta_list = []
    for i in range(0, q1_list.shape[0]):
        quat1 = [q1_list[i], q2_list[i], q3_list[i]]
        mtx_list.append(testroll.quat2mtx(quat1))

    # psai = 0
    # theta = 0
    # phi = 0
    # rot_quat = quat.Quat((psai, theta, phi))
    x, y, z = testroll.sph2cart(ra * np.pi / 180, dec * np.pi / 180, 1)
    xyz = [x, y, z]
    # dcm_b_f_rot = rot_quat.transform.T

    print('length: lighr curvefile, Quat file >>>', len(lctime), len(q1_list))

    atthd.close()
    lchd.close()
    del atthd, lchd


    mtx_list = []
    accept_dex = []
    accept_time = []
    accept_alpha = []
    accept_beta = []
    roll = box_roll[select_box_dex]
    dcm_bx_to_b1 = dcm_roll_theta(roll)
    for i in range(0,q1_list.shape[0]):
        quat1 = [q1_list[i],q2_list[i],q3_list[i]]
        mtx = testroll.quat2mtx(quat1)
        mtx_list.append(mtx)
        delta_alfa0, delta_beta0, xr, yr, zr = testroll.quat2box(dcm_bx_to_b1,dcm_b_f, mtx, xyz)

        if delta_alfa0 < alpha_uplim[select_box_dex] and delta_alfa0 > alpha_downlim[select_box_dex] and delta_beta0 < beta_uplim[select_box_dex] and delta_beta0 > beta_downlim[select_box_dex]:
            flag = 1
            #tot_flag = np.logical_and(np.logical_and(flag[0],flag[1]),flag[2])
            #if tot_flag:
            accept_dex.append(i)
            accept_time.append(lctime[i])
            accept_alpha.append(delta_alfa0)
            accept_beta.append(delta_beta0)

    #-------------------用以多段好时间的筛选-----------------
    print(accept_time,len(accept_time))
    merge_scal = 5
    tem_accept_time = np.append([0],accept_time[:-1])
    print(tem_accept_time,len(tem_accept_time))
    gti_condi = np.greater(accept_time-tem_accept_time, merge_scal)
    print(gti_condi,len(gti_condi),np.sum(gti_condi))
    accept_time = np.array(accept_time)
    gti_time_up = accept_time[gti_condi]
    if len(gti_time_up)>1:
        gti_time_down = np.append(gti_time_up[1:],accept_time[-1])
    elif len(gti_time_up)==1:
        gti_time_down = np.array([accept_time[-1]])
    else:
        print("No accept god time interval!!")

    print(gti_time_up,'\n',gti_time_down)


    boxfile = pf.open('%s' % (evtfile_list[select_box_dex]))
    t = boxfile[1].data.field(0)
    cts = boxfile[1].data.field(1)
    err = boxfile[1].data.field(2)
    ###bkg = boxfile[1].data.field('bkg')
    mask = np.in1d(t, accept_time)
    print("before select:",len(t))
    print("after select:",len(accept_time))
    col1 = pf.Column(name='Time', format='D', array=t[mask])
    col2 = pf.Column(name='counts', format='D', array=cts[mask])
    col3 = pf.Column(name='error', format='D', array=err[mask])
    ###col4 = pf.Column(name='bkg', format='D', array=bkg[mask])
    cols1 = pf.ColDefs([col1, col2, col3])###pf.ColDefs([col1, col2, col3, col4])
    col1 = pf.Column(name='TSTART', format='D', array=gti_time_up)
    col2 = pf.Column(name='TSTOP', format='D', array=gti_time_down)
    cols2 = pf.ColDefs([col1, col2])
    tbhdu1 = pf.BinTableHDU.from_columns(cols1, name='LC')
    tbhdu2 = pf.BinTableHDU.from_columns(cols2, name='GTI')
    prihdu = pf.PrimaryHDU()
    tbhdulist =  pf.HDUList([prihdu, tbhdu1, tbhdu2])###pf.HDUList([prihdu, tbhdu1, tbhdu2])
    lc_file = '%s/new%s_%s_box%s.fits' % (dir_front,instr, year, select_box_dex)
    tbhdulist.writeto(lc_file, clobber=True)
    lc_list.append(lc_file)

    np.save('%s_alpha_box%s' % (instr,str(select_box_dex)), accept_alpha)
    np.save('%s_beta_box%s' % (instr,str(select_box_dex)), accept_beta)
    np.save('%s_cts_box%s' % (instr, str(select_box_dex)), cts[mask])
    np.save('%s_err_box%s' % (instr, str(select_box_dex)), err[mask])
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 1, 1)  # , projection='3d')
    plt.plot(np.array(accept_alpha), np.array(accept_beta), "*")
    fig.savefig('%s_data_position_box%s.png' % (instr,str(select_box_dex)))

# att_list = [infile]
# gen_config(att_list,lc_list,dir_front)


