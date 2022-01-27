import numpy as np

def dex1todex2(n,arr1,arr2):
    tem_dex =np.where(arr1==n)[0]#len
    return arr2[tem_dex]

lefile = 'burst_P021406400301le_1_6.txt'
mefile = 'burst_P021406400301me_7_12.txt'
tot_time_file = 'burst_P021406400301_all.txt'
time_start_flie = 'time_start.txt'
letime_list = np.loadtxt(lefile,delimiter=',')
metime_list = np.loadtxt(mefile,delimiter=',')
time_start = np.loadtxt(time_start_flie)
letime_add_dex = letime_list > -1
metime_add_dex = metime_list > -1
letime_list[letime_add_dex] += time_start[0]
metime_list[metime_add_dex] += time_start[-1]
##print(letime_list[:,0])
#print(metime_list)

meburst = metime_list[:,2:4]
leburst = letime_list[:,2:4]
print(leburst)
print(meburst)

num_le = leburst.shape[0]
num_me = meburst.shape[0]
start_tem_dex = np.zeros((num_le,num_me))
#start_tem_dex2 = np.zeros((num_le,num_me))
stop_tem_dex = np.zeros((num_le,num_me))
#stop_tem_dex2 = np.zeros((num_le,num_me))
start_me_dexplus = np.zeros(num_le)

bin_tem_dex = np.zeros((num_le,num_me))
for i in range(num_le):
    start_tem_dex[i,:] = np.logical_and((meburst[:,0] <= leburst[i,0]), (meburst[:,1] >= leburst[i,0]))
    stop_tem_dex[i,:] = np.logical_and((meburst[:,0] <= leburst[i,1]), (meburst[:,1] >= leburst[i,1]))
    bin_tem_dex[i,:] = np.logical_and((meburst[:,0] >= leburst[i,0]), (meburst[:,1] <= leburst[i,1]))
    condi = np.logical_xor(meburst[:, 0] >= leburst[i, 0], np.append(meburst[1:, 0], 1e10) >= leburst[i, 0])
    condi_sum = np.sum(condi)
    ###print(condi)
    if condi_sum>0:
        start_me_dexplus[i] = np.arange(0,num_me)[condi] + 1


###print(start_me_dexplus)
###-----------------------
# 经过这么操作后，start_tem_dex[i,j]是LE第i个暴开始时间是否在ME第j个暴里，stop_tem_dex[i,j]是LE第i个暴结束时间是否在ME第j个暴里
# start_me_dexplus是LE第i个暴插往哪个ME的后面
###-----------------------


burst_allin = np.logical_and(start_tem_dex, stop_tem_dex)
burst_halfin = np.logical_xor(start_tem_dex, stop_tem_dex)
burst_ifin = np.logical_or(np.logical_or(start_tem_dex, stop_tem_dex),bin_tem_dex)
#print(burst_allin,burst_halfin,burst_ifin)

check_burst_allin = np.sum(burst_allin,axis=1)
check_burst_halfin = np.sum(burst_halfin,axis=1)
check_burst_ifin = np.sum(burst_ifin,axis=1)
###print(check_burst_allin,check_burst_halfin,check_burst_ifin)
'''
num_burst_allin = np.sum(check_burst_allin)
num_burst_halfin = np.sum(check_burst_halfin)
num_burst_ifin = np.sum(check_burst_ifin)
###print(num_burst_allin,num_burst_halfin,num_burst_ifin)
'''
###-----------------------
# burst_ifin[i,j]是说明LE的第i个暴是ME的第j个暴
###-----------------------





check_burst_ifin_bool = (1 - check_burst_ifin).astype(np.bool)
print(check_burst_ifin_bool)
me_dexplus = start_me_dexplus[check_burst_ifin_bool]
print(me_dexplus)


me_insert_dex = np.zeros(num_me)
sub_me_dexplus = np.append(np.append(0,me_dexplus),num_me-1).astype(np.int)
print("LE不在ME中的暴，在ME相应位置对应的指标",me_dexplus)
print(sub_me_dexplus)
for k in range(len(sub_me_dexplus)-1):
    for j in range(sub_me_dexplus[k],sub_me_dexplus[k+1]+1):
            me_insert_dex[j] = j + k
me_insert_dex = me_insert_dex.astype(np.int)


only_le_dex = np.zeros(len(me_dexplus))
for i in range(len(me_dexplus)):
    only_le_dex[i] = int(me_dexplus[i] + i)
only_le_dex = only_le_dex.astype(np.int)
print(only_le_dex)
le_insert_dex = np.zeros(num_le)

me_dex_lein = np.ones(num_le)*999
for i in range(num_le):
    if check_burst_ifin[i]:
        me_dex_lein[i] = np.where(burst_ifin[i,:]==1)[0].astype(np.int)
        le_insert_dex[i] = dex1todex2(me_dex_lein[i],np.arange(num_me),me_insert_dex)

print("LE暴依次在ME中的指标，不在的为999：",me_dex_lein)
check_burst_ifin_bool = (1 - check_burst_ifin).astype(np.bool)
print(check_burst_ifin_bool)
le_dexplus = np.arange(0,num_le)[check_burst_ifin_bool]
print("LE数组里，不在ME的指标",le_dexplus)
le_insert_dex[le_dexplus] = only_le_dex
le_insert_dex = le_insert_dex.astype(np.int)
print("ME第i个时间插入最终数组的指标：",me_insert_dex)
print("LE第i个时间插入最终数组的指标：",le_insert_dex)



num_extent = len(le_dexplus)
num_tot = num_me + num_extent
print(num_extent,num_tot)
tot_time_list = np.ones((num_tot,12))*(-1)
for i in range(num_le):
    tot_time_list[le_insert_dex[i],0:6] = letime_list[i,:]
for j in range(num_me):
    tot_time_list[me_insert_dex[j],6:12] = metime_list[j,:]


# dex = np.where(tot_time_list_tem==0)
# #print(dex[0],dex[1])
# tot_time_list = tot_time_list_tem
# tot_time_list[dex[0],dex[1]] = -1
#print(tot_time_list)
np.savetxt(tot_time_file, tot_time_list, delimiter=',')



# for i in range(num_le):
#     if np.sum(burst_ifin[i, :]) == 0:###LE burst不在ME里, 需要插入
#         k += 1
#     else:
#         for j in range(num_me):
#             tot_time_list_tem[j + k, 6:12] = metime_list[j, :]
#             if j<start_me_dexplus[i]:
#                 if burst_ifin[i,j] == 1:
#                     tot_time_list_tem[j+k-1, 0:6] = letime_list[i,0:6]
#             else:
#                 tot_time_list_tem[j + k, 0:6] = metime_list[j, :]
#                 if burst_ifin[i,j] == 1:
#                     tot_time_list_tem[j+k, 0:6] = letime_list[i, 0:6]