import datetime
import time
import os
import glob
import numpy as np
from multiprocessing import Process
import sys
import smtplib
from email.mime.text import MIMEText
sys.path.append("/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc")
###sys.Wpath.append("/sharefs/hbkg/user/saina/tar/P0201013")
import timing_run
import write_007xml_index

def mkdir_try(dirname):
    if os.path.exists(dirname) == 0:
        try:
            os.makedirs(dirname)
        except OSError:
            print('Wrong with make dir:\n'+dirname)

def file_exists(filename):
    flag = int(os.path.exists(filename))
    return flag

def exec_copy(Wblindlc_2,oWblindlc_2):
    Wblindlc_2_flag = file_exists(Wblindlc_2)
    oWblindlc_2_flag = file_exists(oWblindlc_2)

    print('result filename: ',Wblindlc_2)
    print('original filename: ',oWblindlc_2)
    if Wblindlc_2_flag == 0:
        print('before cp no result file', Wblindlc_2_flag)
        os.system('cp %s %s' % (oWblindlc_2, Wblindlc_2))
    if oWblindlc_2_flag == 0:
        print('No original file!', oWblindlc_2_flag)

    rWblindlc_2_flag = file_exists(Wblindlc_2)

    if rWblindlc_2_flag == 0:
        print('Wrong with execute copy!', rWblindlc_2_flag)
    else:
        print('Copy Done!', rWblindlc_2_flag)

def cpfile(Wpath,ObsID,oWpath = '/sharefs/hbkg/user/saina/data294/P0211007/'):
    OrgWpath = Wpath + '/Org'  ###
    NetWpath = Wpath + '/Net'  ###
    MiddWpath = Wpath + '/Midd'
    ProdWpath = Wpath + '/Prod'
    mkdir_try(OrgWpath)
    mkdir_try(NetWpath)
    mkdir_try(MiddWpath)
    mkdir_try(ProdWpath)
    Org_obspath = "%s/%s" % (OrgWpath, ObsID)  ###
    Net_obspath = "%s/%s" % (NetWpath, ObsID)  ###
    Org_attpath = "%s/%s" % (OrgWpath, 'Att')  ###
    Org_gtipath = "%s/%s" % (OrgWpath, 'GTI')  ###
    Org_ehkpath = "%s/%s" % (OrgWpath, 'EHK')  ###
    Org_pipath = "%s/%s" % (OrgWpath, 'PI')  ###
    Org_screenpath = "%s/%s" % (OrgWpath, 'Screen')  ###
    Midd_obspath = Wpath + '/Midd/' + ObsID  ###
    Prod_obspath = Wpath + '/Prod/' + ObsID
    mkdir_try(Org_obspath)
    mkdir_try(Net_obspath)
    mkdir_try(Org_attpath)
    mkdir_try(Org_gtipath)
    mkdir_try(Org_ehkpath)
    mkdir_try(Org_pipath)
    mkdir_try(Org_screenpath)
    mkdir_try(Midd_obspath)
    mkdir_try(Prod_obspath)

    ###--------------------mkdirs---------------------
    ###--------------------outfiles-------------------
    Wprodfits = "%s/%s_FITDATA_ME.fits"% (Prod_obspath, ObsID)
    Wprodsrc = "%s/src_%s.txt"% (Prod_obspath, ObsID)
    Wattfile = "%s/%s_Att.fits" % (Org_attpath, ObsID)  ###perfect
    Wgtifile = "%s/%s_gtiv2.fits" % (Org_gtipath, ObsID)  ###
    centre_txtfile = "%s/%s/centre.dat" % (OrgWpath, ObsID)  ###
    Wpifile = "%s/%s_me_pi.fits" % (Org_pipath, ObsID)  ###perfect
    Wdeadfile = "%s/%s_me_dead.fits" % (Org_obspath, ObsID)  ###dead time file
    Wgradefile = "%s/%s_me_grade.fits" % (Org_obspath, ObsID)  ###
    Wehkfile = "%s/%s_ehk.fits" % (Org_ehkpath, ObsID)  ###perfect
    Wme_gtifile = "%s/%s_me_gti.fits" % (Org_gtipath, ObsID)  ###perfect
    Wbaddetfile = "%s/%s_me_bdet.fits" % (Org_obspath, ObsID)  ###bad det
    Wscreenfile = "%s/%s_me_screen.fits" % (Org_screenpath, ObsID)  ###
    Wblindlc_2 = "%s/%s/me_2_blind_g0_46.lc" % (NetWpath, ObsID)
    Wsmalllc_2 = "%s/%s/me_lc_box2_small_cut.fits"% (NetWpath, ObsID)
    ###--------------------outfiles-------------------

    oWprodfits = "%s/%s/ME/result/%s_FITDATA.fits"%(oWpath,ObsID,ObsID)
    oWprodsrc = "/hxmt/work/HXMT_scan_data/src_list/ME/src_%s.txt"%(ObsID)
    oWattfile = "%s/%s/att.fits"%(oWpath,ObsID)
    oWgtifile = "%s/%s/ME/gtiv2.fits"%(oWpath,ObsID)
    ocentre_txtfile = "%s/%s/ME/centre.dat"%(oWpath,ObsID)
    oWpifile = "%s/%s/ME/me_pi.fits"%(oWpath,ObsID)
    oWdeadfile = "%s/%s/ME/me_dt.fits"%(oWpath,ObsID)
    oWgradefile = "%s/%s/ME/me_grade.fits"%(oWpath,ObsID)
    oWehkfile = "%s/%s/ehk.fits"%(oWpath,ObsID)
    oWme_gtifile = "%s/%s/ME/gtiv2.fits"%(oWpath,ObsID)
    oWbaddetfile = "%s/%s/bd_%s.fits"%(oWpath,ObsID,ObsID)
    oWscreenfile = "%s/%s/ME/me_screen.fits"%(oWpath,ObsID)
    oWblindlc_2 = "%s/%s/ME/me_2_blind.lc"%(oWpath,ObsID)
    oWsmalllc_2 = "%s/%s/ME/me_lc_box2_small_cut.fits"%(oWpath,ObsID)
    ###
    ###exec_copy(Wblindlc_2,oWblindlc_2)
    ###exec_copy(Wsmalllc_2,oWsmalllc_2)
    ###exec_copy(Wbaddetfile,oWbaddetfile)
    exec_copy(Wprodfits,oWprodfits)
    exec_copy(Wprodsrc,oWprodsrc)

def main(h=0, m=0):
    while True:
        stime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
        # obin='P0101294'
        if len(sys.argv) < 2:
            print
            "Need the ffn!"
            obin = 'P0' + str(input("ff ObsID (211007 211009 201013 101294 101295): "))
            print
            obin
        # sys.exit(1)
        else:
            print
            sys.argv[1]
            obin = sys.argv[1]
        fpath = "/hxmt/work/HXMT-DATA/1L/A%s/%s" % (obin[1:3], obin)
        '''
        obpath='/hxmt/work/HXMT-DATA/1L/'
        if obin[-1]=='7':
            fpath=obpath+'A02/'+obin###+'/'
        if obin[-1]=='4':
            fpath=obpath+'A01/'+obin###+'/'
        if obin[-1]=='9':
            fpath=obpath+'A02/'+obin###+'/'
        if obin[-1]=='5':
            fpath=obpath+'A01/'+obin###+'/'
        if obin[-1]=='3':
            fpath=obpath+'A02/'+obin###+'/'
        '''
        ldata_olist = glob.glob(r'%s/%s*/%s*/' % (fpath, obin, obin))
        ldatalist = []
        for i in ldata_olist:
            ldatalist = np.append(ldatalist, i[-29:-16])

        Wpath = '/sharefs/hbkg/data/SCAN/ME'
        if obin == 'P0101294' or obin == 'P0201013' or obin == 'P0211007' or obin == 'P0211009' or obin == 'P0301240' or obin == 'P0311313' or obin == 'P0311314':
            oWpath = '/sharefs/hbkg/user/saina/data294' + '/%s' % obin
        elif obin == 'P0101295':
            oWpath = '/sharefs/hbkg/user/saina/data295' + '/%s' % obin
        ###----------------------------------position for genlc------------------------------------
        '''
        if obin[-1]=='5':
            lcpath='/sharefs/hbkg/user/saina/data295/'+obin+'/'
        else:
            lcpath='/sharefs/hbkg/user/saina/data294/'+obin+'/'
        '''
        scan_path = '/sharefs/hbkg/data/SCAN/ME'
        lcpath = '/sharefs/hbkg/data/SCAN/ME/Net'

        lclist = os.listdir(lcpath)
        # ldatalist=os.listdir("/hxmt/work/HXMT-DATA/1L/A02/P0211007/")
        lclistre = []
        for i in lclist:
            if not os.path.exists(lcpath + "/%s/me_lc_box2_small_cut.fits" % (str(i))):
                lclistre = np.append(lclistre, i)
        prodlist = []
        for j in ldatalist:
            if os.path.exists(scan_path + "/Prod/%s/%s_FITDATA_ME.fits" % (str(j), str(j))):
                prodlist = np.append(prodlist, j)

        # ldatalist=ldatalist[:700]
        # print("ldatalist",ldatalist)
        # mark=np.in1d(ldatalist,lclist,invert=True)
        # print("lclist",lclist)
        # print("mark",mark)
        # ldatalist=np.array(ldatalist)
        mins = list(set(ldatalist) - set(prodlist))
        ndatalist = np.r_[mins, lclistre]
        ndatalist.sort()
        print("ndatalist", ndatalist)
        list5 = [[]] * int(1 + len(ndatalist) / 5)
        for i in range(len(list5)):
            list5[i] = ndatalist[i * 5:(i + 1) * 5]
        # print list5[i]

        for i in ndatalist:
            cpfile(Wpath, str(i), oWpath)


main(0, 0)