import datetime
import time
import os
import glob
import numpy as np
from multiprocessing import Process
import sys
import smtplib
from email.mime.text import MIMEText

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
    mkdir_try(OrgWpath)
    mkdir_try(NetWpath)
    mkdir_try(MiddWpath)
    Org_obspath = "%s/%s" % (OrgWpath, ObsID)  ###
    Net_obspath = "%s/%s" % (NetWpath, ObsID)  ###
    Org_attpath = "%s/%s" % (OrgWpath, 'Att')  ###
    Org_gtipath = "%s/%s" % (OrgWpath, 'GTI')  ###
    Org_ehkpath = "%s/%s" % (OrgWpath, 'EHK')  ###
    Org_pipath = "%s/%s" % (OrgWpath, 'PI')  ###
    Org_screenpath = "%s/%s" % (OrgWpath, 'Screen')  ###
    Midd_obspath = Wpath + '/Midd/' + ObsID  ###
    mkdir_try(Org_obspath)
    mkdir_try(Net_obspath)
    mkdir_try(Org_attpath)
    mkdir_try(Org_gtipath)
    mkdir_try(Org_ehkpath)
    mkdir_try(Org_pipath)
    mkdir_try(Org_screenpath)
    mkdir_try(Midd_obspath)
    ###--------------------mkdirs---------------------
    ###--------------------outfiles-------------------
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
    exec_copy(Wsmalllc_2,oWsmalllc_2)
    exec_copy(Wbaddetfile,oWbaddetfile)


def main(h=0, m=0):
    while True:
        stime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
        # obin='P0101294'
        if len(sys.argv) < 2:
            print("Need the ffn!")
            obin = 'P0' + str(input("ff ObsID (211007 211009 201013 101294 101295): "))
            print(obin)
        # sys.exit(1)
        else:
            print(sys.argv[1])
            obin = sys.argv[1]

        fpath = "/hxmt/work/HXMT-DATA/1L/A%s/%s" % (obin[1:3], obin)
        '''
        obpath = '/hxmt/work/HXMT-DATA/1L/'
        if obin[-1] == '7':
            fpath = obpath + 'A02/' + obin  ###+'/'
        if obin[-1] == '4':
            fpath = obpath + 'A01/' + obin  ###+'/'
        if obin[-1] == '9':
            fpath = obpath + 'A02/' + obin  ###+'/'
        if obin[-1] == '5':
            fpath = obpath + 'A01/' + obin  ###+'/'
        if obin[-1] == '3':
            fpath = obpath + 'A02/' + obin  ###+'/'
        '''
        ldata_olist = glob.glob(r'%s/%s*/%s*/' % (fpath, obin, obin))
        ldatalist = []
        for i in ldata_olist:
            ldatalist = np.append(ldatalist, i[-29:-16])
        print('ldatalist:',ldatalist)

        Wpath = '/sharefs/hbkg/data/SCAN/ME'
        if obin== 'P0101294' or obin== 'P0201013' or obin=='P0211007' or obin=='P0211009' or obin=='P0301240' or obin=='P0311313' or obin=='P0311314':
            oWpath = '/sharefs/hbkg/user/saina/data294'
        elif obin=='P0101295':
            oWpath = '/sharefs/hbkg/user/saina/data295'
        for i in ldatalist:
            cpfile(Wpath, str(i), oWpath+'/%s'%obin)
        ###----------------------------------position for genlc------------------------------------
        '''
        if obin[-1]=='5':
            lcpath='/sharefs/hbkg/user/saina/data295/'+obin+'/'
        else:
            lcpath='/sharefs/hbkg/user/saina/data294/'+obin+'/'
        '''
        lcpath = '/sharefs/hbkg/data/SCAN/ME/Net'

        lclist = os.listdir(lcpath)
        # ldatalist=os.listdir("/hxmt/work/HXMT-DATA/1L/A02/P0211007/")
        lclistre = []

        for i in lclist:
            #cpfile(Wpath, str(i), oWpath+'/%s'%obin)
            if not os.path.exists(lcpath + "/%s/me_lc_box2_small_cut.fits" % (str(i))):
                lclistre = np.append(lclistre, i)

        # ldatalist=ldatalist[:700]

        mark = np.in1d(ldatalist, lclist, invert=True)
        ldatalist = np.array(ldatalist)
        ndatalist = np.r_[ldatalist[mark], lclistre]
        ndatalist.sort()
        print(ndatalist,len(ndatalist))
        while True:
            now = datetime.datetime.now()
            # if (now.hour%3)==0 and (now.minute%5)==m:
            if (now.minute % 5) == m:
                break
            time.sleep(60)
        '''
        list5 = [[]] * int(1 + len(ndatalist) / 5)
        for i in range(len(list5)):
            list5[i] = ndatalist[i * 5:(i + 1) * 5]
        # print list5[i]

        rlist = []
        doSth(list5, fpath, lcpath)
        for k in list5:
            for num in k:
                if os.path.exists(lcpath + "/%s/me_lc_box2_small_cut.fits" % (str(num))):
                    rlist = np.append(rlist, num)

        # rlist = os.listdir("/sharefs/hbkg/user/saina/data294/P0211007/")
        # mark2=np.in1d(rlist,lclist,invert=True)
        # aplist=[]
        # if True in mark2:
        # aplist=np.array(rlist)[mark2]
        aplist = rlist
        # stime=time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
        title = 'New data list: %s' % stime
        if len(aplist) > 0:
            content = '%s' % aplist
        #	send_mail(email_user, email_pwd, maillist, title, content)
        else:
            content = 'No new data today'
        ###change apfile by LQ
        ###apfile = '/sharefs/hbkg/user/saina/tar/'+obin+'/aplist.txt'
        apfile = '/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/apfiles/' + obin + '/aplist.txt'
        with open(apfile, 'a') as f:
            for i in aplist:
                f.write('%s' % i + '\n')
        '''

main(0,0)