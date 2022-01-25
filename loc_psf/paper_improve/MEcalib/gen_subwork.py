listall = []
dir_strlist = []
with open('dir_config.txt', 'r') as f:
    for line in f.readlines():
        line = line.strip()
        dir_strlist.append(line)
print(dir_strlist)
program_tree = dir_strlist[0]
scan_tree = dir_strlist[1]

with open('crablist.txt','r')as f:
    for i in f:
        listall.append(i[:13])

for i in range(len(listall)):
    ObsID = listall[i]
    command1 = 'sh sub_gps_genlc.sh %s %s'%(ObsID,'7_12')
    command2 = 'sh sub_gps_genlc.sh %s %s'%(ObsID,'7_20')
    with open(program_tree + '/calib7_12_mission.sh', 'a+') as f:
        print >> f, command1
    with open(program_tree + '/calib7_20_mission.sh', 'a+') as f:
        print >> f, command2