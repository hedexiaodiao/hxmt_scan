import sys,os
import shutil

def mkdir_try(dirname):
    if os.path.exists(dirname) == 0:
        try:
            os.makedirs(dirname)
        except OSError:
            print('Wrong with make dir:\n'+dirname)

def all_dir(mkdir_path,scan_dir):
    pyfiles_dir = mkdir_path + '/HXMT_GPS_Program'
    mkdir_try(pyfiles_dir)

    HE_gps_sml = scan_dir + '/HE_GPS_SML'
    HE_GRB = scan_dir + '/HE_GRB'
    HE_data = scan_dir + '/HE_data'
    HE_Production = scan_dir + "/HE_Production"
    HE_srclist = HE_Production + '/srclist'
    mkdir_try(scan_dir)
    mkdir_try(HE_gps_sml)
    mkdir_try(HE_GRB)
    mkdir_try(HE_data)
    mkdir_try(HE_Production)
    mkdir_try(HE_srclist)

    PSF = scan_dir + '/PSF'
    subpsf = PSF + '/psf'
    mkdir_try(PSF)
    mkdir_try(subpsf)
    shutil.copyfile('/sharefs/hbkg/user/luoqi/GRB/work/ihep4/HXMT_GPS_Program/config_he.xml', subpsf+'/config_he.xml')

    ort = scan_dir + '/OrtJ2000'
    mkdir_try(ort)

    pfile = mkdir_path + '/pfiles'
    mkdir_try(pfile)

    datapath = HE_data

    sublcpath = datapath + '/lc_1811/'
    grblcpath = datapath + '/lc_grb/'
    mkdir_try(sublcpath)
    mkdir_try(grblcpath)

    ehkpath = datapath + '/ehk/'
    pipath = datapath + '/pi/'
    gtipath= datapath + '/gti/'
    screenpath = datapath + '/screen/'
    crtattpath = datapath + '/att/'
    mkdir_try(ehkpath)
    mkdir_try(pipath)
    mkdir_try(gtipath)
    mkdir_try(screenpath)
    mkdir_try(crtattpath)

    lcpath = sublcpath + '/lc/'
    netlcpath=  sublcpath + '/netlc/'
    mkdir_try(lcpath)
    mkdir_try(netlcpath)

if __name__=="__main__":
    dir_strlist = []
    with open('dir_config.txt', 'r', encoding='utf-8') as f:
        for line in f.readlines():
            line = line.strip()
            dir_strlist.append(line)
    print(dir_strlist)
    mkdir_path = dir_strlist[0]
    scan_dir = dir_strlist[1]
    if len(sys.argv) < 2:
        pass
        ###mkdir_path = '..'
        ###scan_dir = '/sharefs/hbkg/data/SCAN/luoqi'
    elif len(sys.argv) < 3:
        mkdir_path = sys.argv[1]
    else:
        mkdir_path = sys.argv[1]
        scan_dir = sys.argv[2]
    all_dir(mkdir_path,scan_dir)