import os
import numpy as np
import glob
'''
for i in range(2,22+1):
    command = 'python genhe.py P0101299008{:02d} note.xml'.format(i)
    os.system(command)
'''
id_end2 = [1,2,3,4,6,7,9,10,11,12,13]
for i in id_end2:
    for j in range(24):
        obsid = 'P0101299{:03d}{:02d}'.format(i,j)
        basedatapath = '/hxmt/work/HXMT-DATA/1L/'+"A%s/%s/%s/"%(obsid[1:3],obsid[:8],obsid[:11])
        ###if os.path.exists(basedatapath):
        ###    print(obsid,'exist!!!!')
        datapath_list = glob.glob(basedatapath+obsid+'-*')
        if len(datapath_list)>0:
            datapath = datapath_list[0]
            if os.path.exists(datapath):
                print(obsid,'exist!!!!')
                command = 'python hegen_v1.py '+obsid
                print(command)
                os.system(command)
