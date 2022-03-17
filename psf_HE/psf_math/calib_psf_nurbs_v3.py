#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: calib_psf_nurbs.py
# @Author: luoqi
# @Institution: Institute of High Energy Physics, CAS
# @E-mail: woshiluoqi123@outlook.com, luoqi@ihep.ac.cn
# @Site:
# @Time: 3月 07, 2022
# ---

from geomdl import NURBS
from geomdl import utilities as utils
from geomdl import compatibility as compat
from geomdl.visualization import VisMPL
import numpy as np
from geomdl import operations

#
# Surface exported from your CAD software
#
#-----------平面距离----------
def dis_surf(a0,b0,a1,b1):
    return ((a0-a1)**2 + (b0-b1)**2)**0.5


def get_pn_vw(x,y,data_x,data_y,data_z,data_zerr,x_bound,n):

    #---------------相当于线性插值-----------
    if n==1:
        # ----------左上--------
        save_dex = np.logical_and(x >= data_x, y <= data_y)
    elif n==2:
        # ----------右上--------
        save_dex = np.logical_and(x <= data_x, y <= data_y)
    elif n==3:
        # ----------右下--------
        save_dex = np.logical_and(x <= data_x, y >= data_y)
    elif n==4:
        # ----------左下--------
        save_dex = np.logical_and(x >= data_x, y >= data_y)
    else:
        print('error with point setting n!')
    save_x = data_x[save_dex]
    save_y = data_y[save_dex]
    save_z = data_z[save_dex]
    save_zerr = data_zerr[save_dex]
    sort_dex = np.argsort(dis_surf(x, y, save_x, save_y))
    if len(save_x)>=2:
        d0 = dis_surf(x, y, save_x[sort_dex[0]], save_y[sort_dex[0]])
        d1 = dis_surf(save_x[sort_dex[0]], save_y[sort_dex[0]], save_x[sort_dex[1]], save_y[sort_dex[1]])
        z1 = save_z[sort_dex[1]]
        z0 = save_z[sort_dex[0]]
        zerr1 = save_zerr[sort_dex[1]]
        zerr0 = save_zerr[sort_dex[0]]
        value = z0 + 1*(z0-z1)
        zerr = (zerr1**2 + 1*(zerr0**2 - zerr1**2))**0.5
        weight = f_tan(d0,x_bound)
    elif len(save_x) ==1:
        d0 = dis_surf(x, y, save_x[sort_dex[0]], save_y[sort_dex[0]])
        z0 = save_z[sort_dex[0]]
        zerr0 = save_zerr[sort_dex[0]]
        value = z0
        zerr = zerr0
        weight = f_tan(d0, x_bound)
    else:
        zerr = 0
        value = 0
        weight = 0
    return value, zerr, weight


#-----------平方反比---------
def inv_sq(x):
    return 1/x**2

def get_bound(instru):
    if instru == 'HE':
        Talfa_bound = 5.7
        Tbeta_bound = 1.1
        dcm_b_f = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    elif instru == 'ME':
        Talfa_bound = 4
        Tbeta_bound = 1
        dcm_b_f = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
    elif instru == 'LE':
        Talfa_bound = 6
        Tbeta_bound = 1.6
        dcm_b_f = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    return Talfa_bound, Tbeta_bound

def f_tan(x,x_bound):
    y = 1 - np.tan(np.deg2rad(x))/np.tan(np.deg2rad(x_bound))
    return y

#-------------长方形沿长短边等距划分-----------
def get_rectangle_divd(a_bound,b_bound,p_size_u,p_size_v):
    delta_x = 2*a_bound/(p_size_u-1)
    delta_y = 2*b_bound/(p_size_v-1)
    ctrl_x = np.zeros(p_size_u)
    ctrl_y = np.zeros(p_size_v)
    for i in range(p_size_u):
        ctrl_x[i] = -a_bound + i*delta_x
    for j in range(p_size_v):
        ctrl_y[j] = -b_bound + j*delta_y
    return ctrl_x, ctrl_y

#---------------给一个点和一列点组，找到最近的四个dex和距离--------
def get_nearst4pts(x0,b0,data_x,data_y):
    dis_all = dis_surf(x0,b0,data_x,data_y)
    dex_array = np.argsort(dis_all)[0:4]
    dis_array = dis_all[dex_array]
    return dis_array,dex_array

def get_z_zerr(x0,y0,data_x,data_y,data_z,data_zerr,a_bound,b_bound):
    x_bound = (a_bound**2 + b_bound**2)**0.5
    value_array = np.zeros(4)
    zerr_array = np.zeros(4)
    weight_array = np.zeros(4)
    for i in range(4):
        value_array[i],zerr_array[i],weight_array[i] =get_pn_vw(x0, y0, data_x, data_y, data_z,data_zerr, x_bound, i+1)

    ##print(value_array,weight_array)

    if np.max(weight_array)==1:
        calc_z = value_array[weight_array==1]
        calc_weight = zerr_array[weight_array==1]
    else:
        calc_z = np.sum(value_array*weight_array)/np.sum(weight_array)
        calc_weight = np.sum(zerr_array*weight_array)/np.sum(weight_array)
    return calc_z,calc_weight

#-------------搜索data中最接近ctrl_xy的四个点，误差平方反比作为权重、距离平方反比插值数值--------
def get_ctrlpts(data_x,data_y,data_z,data_zerr,ctrl_x,ctrl_y,a_bound,b_bound):
    lenx = len(ctrl_x)
    leny = len(ctrl_y)
    p_ctrlpts = np.zeros((lenx*leny,3))
    p_weights = np.zeros(lenx*leny)
    for i in range(lenx):
        for j in range(leny):
            ctrl_dex = j*lenx + i
            calc_z, calc_weight = get_z_zerr(ctrl_x[i],ctrl_y[j],data_x,data_y,data_z,data_zerr,a_bound,b_bound)
            p_ctrlpts[ctrl_dex,0] = ctrl_x[i]
            p_ctrlpts[ctrl_dex,1] = ctrl_y[j]
            p_ctrlpts[ctrl_dex,2] = calc_z
            p_weights[ctrl_dex] = calc_weight
    return p_ctrlpts, p_weights

def new_get_ctrlpts(data_x,data_y,data_z,data_zerr,a_bound,b_bound):
    p_ctrlpts = np.zeros((len(data_x),3))
    p_ctrlpts[:,0] = data_x
    p_ctrlpts[:,1] = data_y
    p_ctrlpts[:,2] = data_z
    p_weights = inv_sq(data_zerr)
    return p_ctrlpts, p_weights


###需要有p_size_u >= p_degree_u+1, p_size_v >= p_degree_v+1, 一个是插值点, 一个是基元幂次数
# # Dimensions of the control points grid
# p_size_u = 4
# p_size_v = 3

# Control points in u-row order
###控制点x、y坐标走向尽量按照等距先x变动再y变动的顺序
# p_ctrlpts = [[0, 0, 0], [1, 0, 3], [2, 0, 0], [3, 0, 0],
#              [0, 1, 0.2], [1, 1, 0.5], [2, 1, 0], [3, 1, -1],
#              [0, 2, -0.6], [1, 2, 0], [2, 2, 3], [3, 2, 0]]
#
# # Weights vector
# p_weights = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

# # Degrees
# p_degree_u = 3
# p_degree_v = 2


#--------------------渲染----------------
def psf_render(instru,box_dex,p_ctrlpts,p_weights,p_size_u, p_size_v,p_degree_u,p_degree_v):
    # Combine weights vector with the control points list
    t_ctrlptsw = compat.combine_ctrlpts_weights(p_ctrlpts, p_weights)

    # Since NURBS-Python uses v-row order, we need to convert the exported ones
    n_ctrlptsw = compat.flip_ctrlpts_u(t_ctrlptsw, p_size_u, p_size_v)

    # Since we have no information on knot vectors, let's auto-generate them
    n_knotvector_u = utils.generate_knot_vector(p_degree_u, p_size_u)
    n_knotvector_v = utils.generate_knot_vector(p_degree_v, p_size_v)

    #
    # Import surface to NURBS-Python
    #

    # Create a NURBS surface instance
    surf = NURBS.Surface()
    print('step 1')

    # Fill the surface object
    surf.degree_u = p_degree_u
    surf.degree_v = p_degree_v
    surf.set_ctrlpts(n_ctrlptsw, p_size_u, p_size_v)
    surf.knotvector_u = n_knotvector_u
    surf.knotvector_v = n_knotvector_v
    print('step 2')
    ###决定了按照u/v的细分程度
    # Set evaluation delta
    surf.delta = 0.025

    ###结点细化在psf标定较为有用
    # Refine knot vectors
    operations.refine_knotvector(surf, [3, 0])

    # Set visualization component
    vis_comp = VisMPL.VisSurface()
    surf.vis = vis_comp
    print('step 3')
    surf_points = surf.evalpts
    print(np.array(surf_points).shape)
    # Render the surface
    # surf.render()
    np.save('%s_psf_box%s' % (instru, str(box_dex)),np.array(surf_points))
    np.savetxt('%s_psf_alpha_box%s.txt' % (instru, str(box_dex)), np.array(surf_points)[:,0])
    np.savetxt('%s_psf_beta_box%s.txt' % (instru, str(box_dex)), np.array(surf_points)[:,1])
    np.savetxt('%s_psf_value_box%s.txt' % (instru, str(box_dex)), np.array(surf_points)[:,2])

def read_data(instru,box_dex):
    data_x = np.load('%s_alpha_box%s.npy' % (instru,str(box_dex)))
    data_y = np.load('%s_beta_box%s.npy' % (instru, str(box_dex)))
    data_z = np.load('%s_cts_box%s.npy' % (instru, str(box_dex)))
    data_zerr = np.load('%s_cts_box%s.npy' % (instru, str(box_dex)))
    return data_x,data_y,data_z,data_zerr

def run_main(instru,box_dex):
    data_x,data_y,data_z,data_zerr = read_data(instru,box_dex)
    a_bound,b_bound = get_bound(instru)

    #-------tem for wrong zerr-------
    data_zerr = np.sqrt(np.abs(data_z))
    condi = data_z > 0
    data_z = np.where(condi, data_z, np.zeros_like(data_z))

    center_z,center_zerr = get_z_zerr(0,0,data_x,data_y,data_z,data_zerr,a_bound,b_bound)
    print('center_z,np.max(data_z)',center_z,np.max(data_z))

    #----------standard model---------
    condi = np.logical_and(data_x < a_bound, data_y < b_bound)
    rad_data_x = np.deg2rad(np.abs(data_x))
    rad_data_y = np.deg2rad(np.abs(data_y))
    rad_a_bound = np.deg2rad(a_bound)
    rad_b_bound = np.deg2rad(b_bound)
    factor = np.where(condi, (1.0 - np.tan(rad_data_x) / np.tan(rad_a_bound)) * \
                              (1.0 - np.tan(rad_data_y) / np.tan(rad_b_bound)),
                              np.zeros_like(rad_data_x))
    factor = factor /(np.sqrt(np.tan(rad_data_x)**2 + np.tan(rad_data_y)**2 + 1))
    #stand_z = np.max(data_z)*factor


    data_z = data_z/np.max(data_z)
    data_z = (data_z-factor) + 0.5*(data_z-factor)**2
    data_zerr = np.ones_like(data_zerr/np.max(data_zerr))
    print(data_z)
    print(np.max(data_z),np.min(data_z),np.mean(data_z))
    print(data_zerr)

    p_size_u = 100
    p_size_v = 100
    ctrl_x, ctrl_y = get_rectangle_divd(a_bound,b_bound,p_size_u,p_size_v)
    p_ctrlpts, p_weights = get_ctrlpts(data_x,data_y,data_z,data_zerr,ctrl_x,ctrl_y,a_bound,b_bound)
    print(np.max(p_ctrlpts[:,2]),np.min(p_ctrlpts[:,2]),np.mean(p_ctrlpts[:,2]))
    p_degree_u = 80
    p_degree_v = 80
    psf_render(instru,box_dex,p_ctrlpts,p_weights,p_size_u, p_size_v, p_degree_u, p_degree_v)

    psf = np.load('%s_psf_box%s.npy' % (instru, str(box_dex)))
    print(np.max(psf[:,2]),np.min(psf[:,2]),np.mean(psf[:,2]))
    print(psf)

instru = 'ME'
for box_dex in range(3):
    run_main(instru,box_dex)

# from geomdl import fitting
# from geomdl.visualization import VisMPL as vis
#
# # The NURBS Book Ex9.1
# points = [[4, 0, 0], [4, 4, 0], [4, 8, 3],[4, 10, 2],[0, 0, 0], [0, 4, 0], [0, 8, -3],[0, 10, -2],
#                   [6, 0, 0], [6, 4, -3], [6, 8, 0],[6, 10, 2],[2, 0, 6], [2, 4, 0], [2, 8, 0],[2, 10, 2],
#                   [8, 0, 2], [8, 4, 1], [8, 8, 1],[8, 10, 2]]
# ##degree = 3  # cubic curve
# size_u = 4
# size_v = 4
# degree_u = 2
# degree_v = 2
#
# # Do global curve approximation
# surf = fitting.approximate_surface(points, size_u,size_v,degree_u,degree_v)
#
# # Plot the interpolated curve
# surf.delta = 0.01
# surf.vis = vis.VisSurface()
# surf.render()
# surf_points = np.array(surf.evalpts)
# print(surf_points,surf_points.shape)