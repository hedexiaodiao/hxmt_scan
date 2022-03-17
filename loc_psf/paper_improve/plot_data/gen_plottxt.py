import numpy as np
import os
# def txtname(txthead,num):
#     return txthead+str(num)+'.txt'
# #txtname = 'plotall_10w_LE_2-6_ME_7-20_burst-'
# txthead = 'plotall_LE_2-6_ME_7-20_burst-'
# #txtname = 'plotall_ME_7-20_burst-'

#-------------------count file----------------------#
def count_dir(data_dir):
    tree1 = os.listdir(data_dir)
    i=0
    dir_list = []
    #print(tree1)
    for element in tree1:
        path = os.path.join(data_dir,element)
        #print(path)
        if not os.path.isdir(path):
            i=i+1
            dir_list.append(path)
    print(dir_list)
    print(i)
    return dir_list,i

if __name__ == '__main__':
    #current_dir = r'G:\ihep5\idl_psf_location\moreLE_2_6_100w'
    #current_dir = r'G:\ihep5\idl_psf_location\LE_2_6_10w'
    #current_dir = r'G:\ihep5\idl_psf_location\moreLE_2_6_ME_7_20_10w'
    #current_dir = r'G:\ihep5\idl_psf_location\moreLE_2_6_ME_7_20_100w'
    #current_dir = r'G:\ihep5\idl_psf_location\right_LE_2-6_ME_7-20_10w'
    ###current_dir = r'G:\ihep5\idl_psf_location\right_LE_2-6_10w'
    ##current_dir = r'G:\ihep5\idl_psf_location\right_LE_2-6_ME_7-12_10w'
    #current_dir = r'G:\ihep5\idl_psf_location\v7_LE_2-6_ME_7-12_10w'
    #current_dir = r'G:\ihep5\idl_psf_location\v7_LE_2-6_ME_7-20_10w'
    #current_dir = r'G:\ihep5\idl_psf_locnew\LE_2-6_ME_7-20_10w'
    #current_dir = r'G:\ihep5\idl_psf_locnew\LE_2-6_ME_7-12_10w'
    current_dir = r'G:\ihep5\idl_psf_locnew\LE_2-6_ME_7-40_10w'
    os.chdir(current_dir)
    txtname_list,tot_num = count_dir('plottxt')
    tem_array = np.loadtxt(txtname_list[0])
    print(tem_array.shape[0])
    write_array = np.zeros((tot_num,tem_array.shape[0]))
    print(write_array.shape[0],write_array.shape[1])
    for i in range(tot_num):
        tem_array = np.loadtxt(txtname_list[i])
        write_array[i,:] = tem_array

    np.savetxt('plotall.txt',write_array)
