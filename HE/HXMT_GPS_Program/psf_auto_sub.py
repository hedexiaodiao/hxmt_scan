import os,re



def shortID():
    Crab_ObsID_no01_list = []
    if os.path.exists('psf_sub_command.sh')==1:
        os.system('mv psf_sub_command.sh psf_sub_command.shbak')
    with open('psf_ObsID.txt', 'r', encoding='utf-8') as f:###get short IDs from f
        with open('psf_sub_command.sh', 'a') as r:###prepare write genlc command to r
            with open('{:s}/psf_ObsID.txt'.format(gps_mission_dir),'a') as t:###write short IDs to t
                for ObsID_no01 in f.readlines():
                    Crab_ObsID_no01_list.append(ObsID_no01)
                    f_no01_path = "/hxmt/work/HXMT-DATA/1L/A%s/%s/%s/"%(ObsID_no01[1:3], ObsID_no01[:8], ObsID_no01[:11])
                    dir_all = os.listdir(f_no01_path)
                    print(dir_all)
                    pattern0 = str(ObsID_no01).strip()
                    ObsID_list = []
                    for i in range(len(dir_all)):
                        print(pattern0,dir_all[i])
                        print(dir_all[i])
                        print(type(pattern0))
                        print(re.search(pattern0,str(dir_all[i]), re.M|re.I))
                        print(pattern0 in dir_all[i])
                        if re.search(str(pattern0),str(dir_all[i]), re.M|re.I) != None:
                            ObsID_list.append(dir_all[:13])
                            print('sh {:s}/HXMT_GPS_Program/sub_gps_task.sh {:s}'.format(program_dir,dir_all[i][:13]),file=r)
                            print('{:s}'.format(dir_all[i][:13]),file=t)

def totID():
    with open('psf_held_ObsID.txt', 'r', encoding='utf-8') as f:###get totIDs which held
        with open('psf_held_command.sh', 'a') as r:###write re sub genlc command to r
            for ObsID in f.readlines():
                print('sh /sharefs/hbkg/user/luoqi/GRB/work/ihep4/HXMT_GPS_Program/sub_gps_task.sh {:s}'.format(
                    ObsID), file=r)

def no_complt_ID():
    needID_list = []
    successID_list = []
    if os.path.exists('{:s}/psf_hand_command.sh'.format(gps_mission_dir))==1:
        os.remove('{:s}/psf_hand_command.sh'.format(gps_mission_dir))

    with open('{:s}/psf_ObsID.txt'.format(gps_mission_dir), 'r', encoding='utf-8') as f:
        for needID in f.readlines():
            needID_list.append(str(needID).strip())

    with open('{:s}/genlc_success.txt'.format(gps_mission_dir), 'r', encoding='utf-8') as l:
        for successID in l.readlines():
            successID_list.append(str(successID).strip())

    with open('{:s}/psf_hand_command.sh'.format(gps_mission_dir), 'a') as w:
        for i,needID in enumerate(needID_list):
            print('i,needID',i,needID)
            flag = 0
            for j,successID in enumerate(successID_list):
                print('j,successID',j,successID)
                print(needID, successID)
                print(re.search(str(needID), str(successID), re.M | re.I))
                if re.search(str(needID), str(successID), re.M | re.I) != None:
                    flag = 1
            if flag == 0:
                print('sh {:s}/HXMT_GPS_Program/sub_gps_task.sh {:s}'.format(program_dir, needID), file=w)


dir_strlist = []
with open('dir_config.txt', 'r', encoding='utf-8') as f:
    for line in f.readlines():
        line = line.strip()
        dir_strlist.append(line)
print(dir_strlist)
program_dir = dir_strlist[0]
scan_dir = dir_strlist[1]

gps_mission_dir = program_dir + '/gps_mission'
shortID()
#no_complt_ID()
#totID()
