import numpy as np

def stick_start_stop(a_start,a_stop):
    sort_a = np.argsort(a_start)
    a_start = np.array(a_start)[sort_a]
    a_stop = np.array(a_stop)[sort_a]
    remain_start = a_start
    remain_stop = a_stop
    print('st',remain_start,remain_stop)
    for i in range(len(a_start)):
        stick_dex = a_stop == a_start[i]
        if np.sum(stick_dex)==1:
            del_start = a_start[i]#np.append(del_start,a_start[i])
            del_stop = a_stop[stick_dex]#np.append(del_stop,a_stop[stick_dex])
            start_dex = np.arange(len(a_start))[a_start == del_start]
            remain_start = np.delete(remain_start, start_dex)
            stop_dex = np.arange(len(a_stop))[a_stop == del_stop]
            remain_stop = np.delete(remain_stop,stop_dex)
        if np.sum(stick_dex)>1:
            print('wrong with remain_short!')
    return remain_start,remain_stop


def remain_short(a_start,a_stop):
    sort_a = np.argsort(a_start)
    print(sort_a,a_start,a_stop)
    a_start = np.array(a_start)[sort_a]
    a_stop = np.array(a_stop)[sort_a]
    remain_a_start = np.array([])
    remain_a_stop = np.array([])
    for i in range(len(a_start)):
        same_dex = a_start == a_start[i]
        if np.sum(same_dex)>=1:
            same_stop = a_stop[same_dex]
            remain_a_stop = np.append(remain_a_stop,np.min(same_stop))
            remain_a_start = np.append(remain_a_start,a_start[i])
    remain_a_start = np.unique(remain_a_start)
    remain_a_stop = np.unique(remain_a_stop)
    return remain_a_start,remain_a_stop

def remain_only(remain_start,remain_stop,start_candi,stop_candi):
    # print('check input',remain_start,remain_stop,start_candi,stop_candi)
    start_in_flag = np.logical_and((remain_start <= start_candi), (remain_stop >= start_candi))
    stop_in_flag = np.logical_and((remain_start <= stop_candi), (remain_stop >= stop_candi))

    if start_candi==stop_candi:
        return remain_start,remain_stop
    elif np.sum(start_in_flag)==1 or np.sum(stop_in_flag)==1:
        start_in_dex = np.arange(len(remain_start))[start_in_flag]
        stop_in_dex = np.arange(len(remain_stop))[stop_in_flag]
        # print('in flag;',start_in_flag,stop_in_flag,remain_start[start_in_dex],remain_stop[stop_in_dex],remain_start,remain_stop)
        remain_start[start_in_dex] = start_candi
        remain_stop[stop_in_dex] = stop_candi
        return remain_start,remain_stop
    elif np.sum(start_in_flag)==0 and np.sum(stop_in_flag)==0:
        remain_start = np.append(remain_start,start_candi)
        remain_stop = np.append(remain_stop,stop_candi)
        # print('out;',start_in_flag,stop_in_flag,remain_start,remain_stop)
        return remain_start,remain_stop
    else:
        print('wrong!!! don\'t allow more than one intersection')
        return remain_start,remain_stop


def calc_intersection(a_start,a_stop,b_start,b_stop):
    #----------------假设不存在自己区间交叉---------------------

    sort_a = np.argsort(a_start)
    sort_b = np.argsort(b_start)
    a_start = np.array(a_start)[sort_a]
    a_stop = np.array(a_stop)[sort_a]
    b_start = np.array(b_start)[sort_b]
    b_stop  = np.array(b_stop)[sort_b]

    # a_start,a_stop = remain_short(a_start,a_stop)
    # #a_stop, a_start = remain_short(a_stop,a_start)
    # b_start,b_stop = remain_short(b_start,b_stop)
    # #b_stop, b_start = remain_short(b_stop,b_start)
    # print('remain;',a_start, a_stop)
    # print('remain;',b_start, b_stop)
    # a_start,a_stop = stick_start_stop(a_start,a_stop)
    # b_start,b_stop = stick_start_stop(b_start,b_stop)
    # print('stick;', a_start, a_stop)
    # print('stick;', b_start, b_stop)
    # a_start, a_stop = remain_short(a_start, a_stop)
    # #a_stop, a_start = remain_short(a_stop, a_start)
    # b_start, b_stop = remain_short(b_start, b_stop)
    # #b_stop, b_start = remain_short(b_stop, b_start)
    # print('stick;',a_start,a_stop)
    # print('stick;',b_start,b_stop)

    if len(a_start) != len(a_stop):
        print('section a has different len')
    if len(b_start) != len(b_stop):
        print('section b has different len')
    a_len = len(a_start)
    b_len = len(b_start)

    #----------每个a区间做如下操作，先找half in , tot in 记为flag，成立后执行如下：
    #----------计算总体np.sum(flag)是否大于1，大于1要分别补；
    #----------对每个flag位置，如果存在a_left<=b_left则补录b_left, a_right>=b_right则补录b_right，否则留a
    #----------因为是求交集，执行完以上便是全部-------
    remain_start = np.array([])
    remain_stop = np.array([])

    for i in range(a_len):
        # a_flag_all_contain = np.zeros(b_len)
        # a_flag_left_contain = np.zeros(b_len)
        # a_flag_right_contain = np.zeros(b_len)

        #------------start留大的，stop留小的-------------
        a_flag_all_contain =  np.logical_and(np.logical_and((b_start >= a_start[i]), (b_stop <= a_stop[i])),b_stop  >= a_start[i])#(b,b)
        a_flag_left_contain = np.logical_and(np.logical_and(( b_start >= a_start[i]), (b_stop >= a_stop[i])),b_start <= a_stop[i])#(b,a)
        a_flag_right_contain = np.logical_and(np.logical_and((b_start <= a_start[i]), (b_stop <= a_stop[i])),b_stop  >= a_start[i])#(a,b)
        a_flag_in_contain = np.logical_and((b_start <= a_start[i]), (b_stop >= a_stop[i]))#(a,a)
        tot_all = np.sum(a_flag_all_contain)
        tot_left = np.sum(a_flag_left_contain)
        tot_right = np.sum( a_flag_right_contain)
        tot_in = np.sum(a_flag_in_contain)
        # print(i,tot_all,tot_left,tot_right,tot_in)
        if tot_all>=1:
            remain_all_start = b_start[a_flag_all_contain]
            remain_all_stop = b_stop[a_flag_all_contain]
            # print('remain_all_start', remain_all_start,remain_all_stop, tot_all)
            for j in range(tot_all):
                remain_start,remain_stop = remain_only(remain_start,remain_stop,remain_all_start[j],remain_all_stop[j])
        if tot_left>=1:
            remain_left_start = b_start[a_flag_left_contain]
            remain_left_stop = a_stop[i]
            # print('remain_left_start', remain_left_start, remain_left_stop, tot_left)
            for j in range(tot_left):
                # print('af j',remain_start,remain_stop,remain_left_start[j], remain_left_stop)
                remain_start,remain_stop = remain_only(remain_start,remain_stop,remain_left_start[j], remain_left_stop)
        if tot_right>=1:
            remain_right_start = a_start[i]
            remain_right_stop = b_stop[a_flag_right_contain]
            # print('remain_right_stop',remain_right_stop,tot_right)
            for j in range(tot_right):
                # print(j)
                remain_start, remain_stop = remain_only(remain_start, remain_stop,  remain_right_start,remain_right_stop[j])
        if tot_in>=1:
            remain_in_start = a_start[i]
            remain_in_stop = a_stop[i]
            # print('remain_in',a_start[i],a_stop[i],tot_in)
            for j in range(tot_in):
                remain_start,remain_stop = remain_only(remain_start,remain_stop, remain_in_start, remain_in_stop)
        # print('remain_start,stop,i',remain_start,remain_stop,i)


    # remain_start = np.unique(remain_start)
    # remain_stop = np.unique(remain_stop)
    sort_dex = np.argsort(remain_start)
    remain_start = remain_start[sort_dex]
    remain_stop = remain_stop[sort_dex]
    # remain_start = np.unique(remain_start)
    # remain_stop = np.unique(remain_stop)
    ### print(remain_start, remain_stop)

    # del_start = np.array([])
    # for i in range(len(remain_start)):
    #     stop_dex = np.arange(len(remain_stop))[remain_stop == remain_start[i]]
    #     remain_stop = np.delete(remain_stop,stop_dex)
    #     if np.sum(stop_dex)>=1:
    #         del_start = np.append(del_start,remain_start[i])
    # print(del_start)
    # for i in range(len(del_start)):
    #     start_dex = np.arange(len(remain_start))[remain_start == del_start[i]]
    #     remain_start = np.delete(remain_start,start_dex)

    #remain_start,remain_stop = remain_short(remain_start,remain_stop)
    ### remain_start,remain_stop = stick_start_stop(remain_start,remain_stop)
    ### remain_start, remain_stop = remain_short(remain_start, remain_stop)
    ### remain_stop, remain_start = remain_short(remain_stop, remain_start)
    # print(remain_start,remain_stop)
    return remain_start,remain_stop

if __name__ == "__main__":
    test_a_start = [2,6,11,15,21,23,25]
    test_a_stop  = [5,10,14,20,23,24,30]
    test_b_start = [1,4,6,10,15,27,30]
    test_b_stop  = [3,6,8,13,22,30,32]

    calc_intersection(test_a_start, test_a_stop, test_b_start, test_b_stop)