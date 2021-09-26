#!/hxmt/soft/Develop/anaconda2/bin/python
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits as pf
import sys,os
from xml.etree.ElementTree import ElementTree,Element
import time
from glob import glob as glob
#from CorrectAtt import correct
###sys.path.append("/sharefs/hbkg/user/nangyi/HXMT_GPS_Tools")
sys.path.append(os.path.abspath(os.path.dirname(__file__))+'/HXMT_GPS_py3Tools')
from exec_cmd import *
###sys.path.append("/sharefs/hbkg/user/luoqi/GRB/work/ihep4/HXMT_GPS_Toolspy3")
import Findfile
from Lctools_201903 import *
time1= time.time()


############################################################################################

def genlc(ELV, DYE_ELV, COR, T_SAA, TN_SAA, SUN_ANGLE, MOON_ANGLE, ANG_DIST,ehkfile,lcfile, lcblind, lcpath, netlcpath,alllcfile,nodeadlc,minpi, maxpi,filestr):
    evtfile,tpfile,attfile,dtfile,hvfile,pmfile,ortehkfile = FindLC_basefile(key,fpath,INST='HE',EHKFILE=OrgEhkTerms)

    print(evtfile)
    print(attfile)

    my_env = os.environ
    my_env['PFILES'] = "{:s}".format(pfilepath)
    print(my_env)

    if NEWDATA:
        correct(attfile,crtattfile)###
        ###########################################################################################
        print("Start \t",key)
        #************************************************************
        if not OrgEhkTerms:
            print("Start Ehk!")
            commond = "hxmtehkgen orbfile="+ ortehkfile +" attfile="+ crtattfile +" outfile="+ ehkfile +" saafile="+ saafile +" step_sec=0.25  clobber=yes" #mean_phi=0 mean_theta=0 mean_psi=0 leapfile="+leapfile+" rigidity="+rigidity+" saafile="+saafile
            cmd_tem = 'source {:s}/HXMT_GPS_Program/hxmtsoft_v2py2.sh;export PFILES="{:s}";export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'.format(program_tree,pfilepath)+commond
            exec_code,exec_log = exec_cmd(cmd_tem,my_env)
            print(cmd_tem)
            print('exec_code:',exec_code)
        else:
            ehkfile = ortehkfile

        #*************************************************************
        commond = "hepical evtfile=" + evtfile + " outfile=" + pifile
        cmd_tem = 'source {:s}/HXMT_GPS_Program/hxmtsoft_v2py2.sh;export PFILES="{:s}";export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'.format(program_tree,pfilepath)+commond
        exec_code, exec_log = exec_cmd(cmd_tem, my_env)
        print(cmd_tem)
        print('exec_code:', exec_code)
        print("Pi end ")
        timepot,flag = he_nozero_pot_CSI(pifile,ehkfile)
        if flag==0:
            with open(grblogfile,'a')as f:
                f.write('%s\t%s\n'%(key,str(time.strftime('%Y.%m.%d-%H:%M',time.localtime(time.time()))) ))
            print('>>>>>GRB Mode<<<<<')
            filestr = grbfilestr
            lcpath = grblcpath + 'lc/'
            netlcpath=  grblcpath + 'netlc/'

            lcfile = lcpath + key +"he"
            lcblind = lcpath + key +"blind"
            alllcfile = lcpath + key + "A_he"
            nodeadlc = lcpath + key + "N_dead"
            allnetlcpath = netlcpath + key +"A_he"
            minpi = [25]*18
            maxpi = [220]*18
        timestart,timestop = timepot[0],timepot[1]

        #*************************************************************
        print("Start GTi!")
        commond = "hegtigen hvfile=" + hvfile + " tempfile=" + tpfile + " outfile=" + gtifile + " ehkfile=" + ehkfile +" pmfile=  defaultexpr=NONE expr='COR3>%s&&ELV>%s&&DYE_ELV>%s&&SAA_FLAG==0&&T_SAA>%s&&TN_SAA>%s&&SUN_ANG>%s&&MOON_ANG>%s' clobber=yes"%(COR,ELV,DYE_ELV,T_SAA,TN_SAA,SUN_ANGLE,MOON_ANGLE)
        cmd_tem = 'source {:s}/HXMT_GPS_Program/hxmtsoft_v2py2.sh;export PFILES="{:s}";export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'.format(
            program_tree, pfilepath) + commond
        exec_code, exec_log = exec_cmd(cmd_tem, my_env)
        print(cmd_tem)
        print('exec_code:', exec_code)
        print("GTi end ")
        del ELV,DYE_ELV,COR,T_SAA,TN_SAA,SUN_ANGLE,MOON_ANGLE,ANG_DIST

        #*************************************************************
        commond = "hescreen evtfile="+ pifile+ " gtifile="+ gtifile +" outfile=" + screenfile + ' userdetid="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17"  baddetfile="" eventtype=1 anticoincidence="yes" minpulsewidth=30 maxpulsewidth=70 clobber=yes'
        cmd_tem = 'source {:s}/HXMT_GPS_Program/hxmtsoft_v2py2.sh;export PFILES="{:s}";export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'.format(
            program_tree, pfilepath) + commond
        exec_code, exec_log = exec_cmd(cmd_tem, my_env)
        print(cmd_tem)
        print('exec_code:', exec_code)
        print("Screen end ")
    else:
        timepot,flag = he_nozero_pot_CSI(pifile,ehkfile)
        if flag==0:
            with open(grblogfile,'a')as f:
                f.write('%s\t%s\n'%(key,str(time.strftime('%Y.%m.%d-%H:%M',time.localtime(time.time()))) ))
            print('>>>>>GRB Mode<<<<<')
            filestr = grbfilestr
            lcpath = grblcpath + 'lc/'
            netlcpath=  grblcpath + 'netlc/'
            lcfile = lcpath + key +"he"
            lcblind = lcpath + key +"blind"
            alllcfile = lcpath + key + "A_he"
            nodeadlc = lcpath + key + "N_dead"
            allnetlcpath = netlcpath + key +"A_he"
            minpi = [25]*18
            maxpi = [220]*18
        timestart,timestop = timepot[0],timepot[1]

    #*************************************************************
    startm = []
    stoptm = []
    thd = pf.open(screenfile)
    for i in range(len(thd)-3):
        startm.append(thd[i+2].data.field(0)[0])
        stoptm.append(thd[i+2].data.field(1)[-1])

    thd.close()
    startm = np.rint(max(startm))+0.5
    stoptm = np.rint(min(stoptm))-0.5
    #************************************************************
    if NEWDATA:
        print("Lc start!")
        commond = "helcgen evtfile=" +screenfile +" outfile=" + lcfile+" deadfile=" + dtfile + ' userdetid="0,4,7,14,17;1,3,8,11,12;5,6,10,13,15" binsize=1 eventtype=1 starttime=%s stoptime=%s minPI=18 maxPI=83 clobber=yes'%(startm,stoptm)
        cmd_tem = 'source {:s}/HXMT_GPS_Program/hxmtsoft_v2py2.sh;export PFILES="{:s}";export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'.format(
            program_tree, pfilepath) + commond + ';rm {:s}/helcgen.par'.format(pfile)
        exec_code, exec_log = exec_cmd(cmd_tem, my_env)
        print(cmd_tem)
        print('exec_code:', exec_code)
        ###os.system('source /hxmt/work/USERS/nangyi/envV2.sh ;export PFILES="%s"; '%(pfilepath)+ commond + ';rm %s/helcgen.par'%pfile)###

        commond = "helcgen evtfile=" +screenfile +" outfile=" + lcblind +" deadfile=" + dtfile + ' userdetid="16" binsize=1 eventtype=1 starttime=%s stoptime=%s minPI=18 maxPI=83 clobber=yes'%(startm,stoptm)
        cmd_tem = 'source {:s}/HXMT_GPS_Program/hxmtsoft_v2py2.sh;export PFILES="{:s}";export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'.format(
            program_tree, pfilepath) + commond + ';rm {:s}/helcgen.par'.format(pfile)
        exec_code, exec_log = exec_cmd(cmd_tem, my_env)
        print(cmd_tem)
        print('exec_code:', exec_code)
        ###os.system('source /hxmt/work/USERS/nangyi/envV2.sh ;export PFILES="%s"; '%(pfilepath)+ commond + ';rm %s/helcgen.par'%pfile)###

    print("Lc end ")

    ############################################################################################
    gtilt = [list(timestart),list(timestop)]
    gtihd = pf.open(gtifile)
    gtilt=list(zip(*gtilt))
    gti = np.rec.array(gtilt,formats='float,float',names='TSTART,TSTOP')

    for i in range(len(gtihd)-2):
        gtihd[i+1].data = pf.FITS_rec.from_columns(gti)

    gtihd.writeto(allgtifile,clobber = 'True')
    gtihd.close()
    del gtihd

    commond = "hescreen " + "evtfile="+ pifile+ " gtifile="+ allgtifile +" outfile=" + allscreenfile + ' userdetid="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17"  baddetfile="" eventtype=1 anticoincidence="yes" minpulsewidth=30 maxpulsewidth=70 clobber=yes'
    cmd_tem = 'source {:s}/HXMT_GPS_Program/hxmtsoft_v2py2.sh;export PFILES="{:s}";export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'.format(
        program_tree, pfilepath) + commond
    exec_code, exec_log = exec_cmd(cmd_tem, my_env)
    print(cmd_tem)
    print('exec_code:', exec_code)
    print("Screen end ")
    #*************************************************************
    print("Lc start!")
    for i in range(18):
        commond = "helcgen evtfile=" + allscreenfile +" outfile=" + alllcfile+" deadfile=" + dtfile + ' userdetid="%s" binsize=1 eventtype=1 starttime=%s stoptime=%s minPI=%s maxPI=%s clobber=yes'%(i,timestart[0]+0.5,timestop[-1]-0.5,minpi[i],maxpi[i])
        cmd_tem = 'source {:s}/HXMT_GPS_Program/hxmtsoft_v2py2.sh;export PFILES="{:s}";export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'.format(
            program_tree, pfilepath) + commond
        exec_code, exec_log = exec_cmd(cmd_tem, my_env)
        print(cmd_tem)
        print('exec_code:', exec_code)
        commond = "helcgen evtfile=" + allscreenfile +" outfile=" + nodeadlc +" deadfile=" + dtfile + ' userdetid="%s" binsize=1 eventtype=1 starttime=%s stoptime=%s minPI=%s maxPI=%s deadcorr=no clobber=yes'%(i,timestart[0]+0.5,timestop[-1]-0.5,minpi[i],maxpi[i])
        cmd_tem = 'source {:s}/HXMT_GPS_Program/hxmtsoft_v2py2.sh;export PFILES="{:s}";export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'.format(
            program_tree, pfilepath) + commond
        exec_code, exec_log = exec_cmd(cmd_tem, my_env)
        print(cmd_tem)
        print('exec_code:', exec_code)

    print("Lc end ")

    ############################################################################################
    os.system("rm %s/*"%(pfile))###
    os.system("rm -rf %s"%(pfile))###
    ############################################################################################
    lcnames = [lcfile+"_g0_0-17.lc",lcfile+"_g1_1-12.lc",lcfile+"_g2_5-15.lc"]
    norm_3netlc(lcnames,lcpath)
    togeth_18lc_v2(key+"A_he",lcpath)
    togeth_18lc_v2(key+"N_dead",lcpath)
    bkg_muti_poly(key+"A_he",lcpath,netlcpath)
    atthd = pf.open(crtattfile)
    attdt=atthd[1].data
    atttm = attdt.field(0)
    attra=attdt.field(1)
    attdec=attdt.field(2)

    def dist(ra1,dec1,ra2,dec2):
        ytp = (np.sin(dec1*np.pi/360.0 - dec2*np.pi/360.0))**2 + np.cos(dec1*np.pi/180.0)*np.cos(dec2*np.pi/180.0)*(np.sin(ra1*np.pi/360.0 - ra2*np.pi/360.0))**2
        ytp = np.arccos(1-2*ytp)*180./np.pi
        return ytp

    attdis =[]
    atttime =[]
    for i,ra1 in enumerate(attra):
        try:
            dec1=attdec[i]
            ra2=attra[i+2]
            dec2=attdec[i+2]
            tp = dist(ra1,dec1,ra2,dec2)
            attdis.append(tp/(atttm[i+2]-atttm[i]))
            atttime.append(atttm[i])
        except:
            pass

    seq = (np.array(attdis)>0.055)&(np.array(attdis)<0.065)
    attimes = atttm[:-2][seq]
    atthd.close()
    del atthd
    for i in range(3):
        if i >1:
            k=5
        else:
            k=i
        lchd = pf.open(glob(lcpath+key+'he_g%s_%s*.lc'%(i,k))[0])
        allnethd = pf.open(netlcpath + key +'A_he_netlc%s.fits'%k)
        nettm = allnethd[1].data.field(0)
        yseq = lchd[1].data.field(2)==1.0
        lctm = lchd[1].data.field(0)[yseq]
        lctm = lctm[np.in1d(lctm,attimes)]
        allnethd[1].data = allnethd[1].data[(np.in1d(nettm,lctm))]
        print('>>>>>>>> timm lenth',allnethd[1].data.shape[0])
        if allnethd[1].data.shape[0]==0:
            with open('Genlc_zero.log','a') as f:
                f.write(key + '\tzero time!\n')
            system.exit(0)
        allnethd.writeto(netlcpath + key + "he_netlc%s.fits"%k,output_verify="ignore",overwrite=True)
        lchd.close()
        allnethd.close()
        del allnethd
        del lchd


    plotnames = [netlcpath+key+"he_netlc0.fits",netlcpath+key+"he_netlc1.fits",netlcpath+key+"he_netlc5.fits"]
    norm_3netlc(plotnames,'hehe')
    plotnetlc(plotnames,netlcpath+key+"he_netlc.png")
    ############################################################################################
    atthd=pf.open(crtattfile)
    attm = atthd[1].data.field(0)
    allnethd = pf.open(netlcpath + key + "he_netlc0.fits")
    lctm = allnethd[1].data.field(0)
    yep = np.in1d(attm,lctm)
    for i in range(1,len(atthd),1):
        atthd[i].data = atthd[i].data[yep]

    try:
        atthd.writeto(attshortfile)
    except :
        atthd.writeto(attshortfile,overwrite=True)
    allnethd.close()
    del allnethd
    atthd.close()
    del atthd
    ############################################################################################
    try:
        os.makedirs(filestr)###mode=0755
    except OSError:
        pass

    tree = ElementTree()
    tree.parse(xmlpath)
    pathnodes = tree.findall("PATH/outpath")
    outnodes = tree.findall("Outfile_Name/outfilename")
    filenodes = tree.findall("Infile_Name/infilename/infilenamelist")
    outfile = filestr + "/config_he.xml"
    print(outfile)
    filenodes[0].text ='\n  ' + netlcpath + key + "he_netlc0.fits \n"
    filenodes[1].text ='\n  ' + netlcpath + key + "he_netlc1.fits \n"
    filenodes[2].text ='\n  ' + netlcpath + key + "he_netlc5.fits \n"
    filenodes[3].text ='\n  ' + crtattfile +' \n'
    pathnodes[0].text ='\n  '+filestr+'  \n'

    tree.write(outfile, encoding="utf-8",xml_declaration=True)

    ###os.system("cp Lcfit.py %s" % (filestr + "/"))
    ###os.system("cp fitdata_plot.py %s" % (filestr + "/"))
    os.chdir(filestr)
    commond = "sh {:s}sub_gps_lcfit.sh {:s}/config_he.xml".format(program_tree + "/HXMT_GPS_Program/",filestr)###"python {:s}Lcfit.py {:s}config_he.xml".format(program_tree + "/HXMT_GPS_Program/",filestr + "/")
    '''###
    cmd_tem = 'source /sharefs/hbkg/user/luoqi/home/mypython;export PFILES="{:s}";export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'.format(
            program_tree, pfilepath) + commond
    exec_code, exec_log = exec_cmd(cmd_tem, my_env)
    print(cmd_tem)
    print('exec_code:', exec_code)
    '''###
    if flag==0:
        commond = ''
    with open(program_tree+'/gps_mission/genlc_success.txt','a+') as f:
        print(key,file=f)
    with open(program_tree+'/gps_mission/lcfit_misson.sh','a+') as f:
        print(commond,file=f)

    print('Fit end!')

    time2 = time.time()
    print("time expended: ",(time2-time1)/3600.0,' hours')
    tmptm = str(time.strftime('%Y.%m.%d-%H:%M',time.localtime(time.time())))
    if flag!=0:
        with open('Genlc_tmp.log','a') as f:
            f.write(key + '\tNoraml LC Generated!  %s\n'%tmptm)


if __name__ == '__main__':
    run = Findfile.findfile()

    dir_strlist = []
    with open('dir_config.txt', 'r', encoding='utf-8') as f:
        for line in f.readlines():
            line = line.strip()
            dir_strlist.append(line)
    print(dir_strlist)
    program_tree = dir_strlist[0]
    scan_tree = dir_strlist[1]

    if len(sys.argv) < 2:
        key = input("input keywords like P010129400101 :")
    elif len(sys.argv)==2:
        key = sys.argv[1]
    else:
        print('Wrong number of parameters!')
        os.exit(0)


    fpath = '/hxmt/work/HXMT-DATA/1L/A%s/%s/' % (key[1:3], key[:-5])
    fpath = "/hxmt/work/HXMT-DATA/1L/A%s/%s/%s/" % (key[1:3], key[:-5], key[:-2])

    ELV, DYE_ELV, COR, T_SAA, TN_SAA, SUN_ANGLE, MOON_ANGLE, ANG_DIST = 6, 6, 0, 5, 5, 30, 7, 360
    minpi = [22, 23, 22, 23, 22, 22, 23, 22, 22, 22, 22, 22, 22, 22, 22, 22, 23, 22]
    maxpi = [80, 85, 81, 85, 83, 84, 87, 82, 82, 83, 82, 82, 82, 81, 80, 82, 83, 82]
    OrgEhkTerms = False
    filestr = scan_tree + "/HE_GPS_SML/" + key[5:]


    grbfilestr = scan_tree + "/HE_GRB/" + key[5:]
    grblogfile = program_tree + "/HXMT_GPS_Program/GRB_DataLog.txt"
    xmlpath = scan_tree + "/PSF/psf/config_he.xml"

    datapath = scan_tree + "/HE_data/"
    sublcpath = datapath + 'lc_1811/'
    grblcpath = datapath + 'lc_grb/'

    ehkpath = datapath + 'ehk/'
    pipath = datapath + 'pi/'
    gtipath = datapath + 'gti/'
    screenpath = datapath + 'screen/'
    crtattpath = datapath + 'att/'

    lcpath = sublcpath + 'lc/'
    netlcpath = sublcpath + 'netlc/'

    crtattfile = crtattpath + key + 'crtAtt.fits'
    attshortfile = crtattpath + key + 'shortAtt.fits'
    ehkfile = ehkpath + key + "ehk.fits"
    pifile = pipath + key + "hepi.fits"
    gtifile = gtipath + key + "hegti.fits"
    screenfile = screenpath + key + "hescreen.fits"

    allgtifile = gtipath + key + "A_hegti.fits"
    allscreenfile = screenpath + key + "A_hescreen.fits"

    lcfile = lcpath + key + "he"
    lcblind = lcpath + key + "blind"
    alllcfile = lcpath + key + "A_he"
    nodeadlc = lcpath + key + "N_dead"
    allnetlcpath = netlcpath + key + "A_he"
    ortpath = scan_tree + '/OrtJ2000/'

    try:
        commond = sys.argv[2]
        NEWDATA = False
    except:
        os.system("rm %s*/" % (datapath) + key + "*")###
        NEWDATA = True

    leapfile = "/hxmt/home/hxmtsoft2/hxmtehkgen/refdata/leapsec.fits"
    rigidity = "/hxmt/home/hxmtsoft2/hxmtehkgen/refdata/rigidity_20060421.fits"
    saafile = "/home/hxmt/guanj/zhaohsV2/hxmtehkgen/SAA/SAA.fits"

    pfile = program_tree + "/pfiles/%s" % str(key)
    print(key)
    '''
    os.system('mkdir %s' % pfile)
    # pfilepath = "%s;/home/hxmt/zhaohaish/soft/install/x86_64-unknown-linux-gnu-libc2.12/syspfiles"%pfile
    pfilepath = "%s;/home/hxmt/hxmtsoft2/soft/install/x86_64-unknown-linux-gnu-libc2.12/syspfiles" % pfile
    os.system("rm %s/*" % pfile)###

    genlc(ELV, DYE_ELV, COR, T_SAA, TN_SAA, SUN_ANGLE, MOON_ANGLE, ANG_DIST,ehkfile, lcfile, lcblind, lcpath, netlcpath, alllcfile,nodeadlc, minpi, maxpi,filestr)
    '''