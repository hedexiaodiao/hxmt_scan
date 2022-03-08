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

def get_z_zerr(x0,b0,data_x,data_y,data_z,data_zerr):
    dis_array,dex_array = get_nearst4pts(x0,b0,data_x,data_y)
    nrst4z= data_z[dex_array]
    nrst4zerr = data_zerr[dex_array]
    if np.min(dis_array)==0:
        calc_z = data_z[dex_array[0]]
        calc_weight = data_zerr[dex_array[0]]
    else:
        calc_z = np.sum(inv_sq(nrst4zerr)*inv_sq(dis_array)*nrst4z)/np.sum(inv_sq(nrst4zerr)*inv_sq(dis_array))
        calc_weight = np.sum(inv_sq(nrst4zerr)*inv_sq(dis_array))/np.sum(inv_sq(dis_array))
    return calc_z,calc_weight

#-------------搜索data中最接近ctrl_xy的四个点，误差平方反比作为权重、距离平方反比插值数值--------
def get_ctrlpts(data_x,data_y,data_z,data_zerr,ctrl_x,ctrl_y):
    lenx = len(ctrl_x)
    leny = len(ctrl_y)
    p_ctrlpts = np.zeros((lenx*leny,3))
    p_weights = np.zeros(lenx*leny)
    for i in range(lenx):
        for j in range(leny):
            ctrl_dex = j*lenx + i
            calc_z, calc_weight = get_z_zerr(ctrl_x[i],ctrl_y[j],data_x,data_y,data_z,data_zerr)
            p_ctrlpts[ctrl_dex,0] = ctrl_x[i]
            p_ctrlpts[ctrl_dex,1] = ctrl_y[j]
            p_ctrlpts[ctrl_dex,2] = calc_z
            p_weights[ctrl_dex] = calc_weight
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

    # Fill the surface object
    surf.degree_u = p_degree_u
    surf.degree_v = p_degree_v
    surf.set_ctrlpts(n_ctrlptsw, p_size_u, p_size_v)
    surf.knotvector_u = n_knotvector_u
    surf.knotvector_v = n_knotvector_v

    ###决定了按照u/v的细分程度
    # Set evaluation delta
    surf.delta = 0.025

    ###结点细化在psf标定较为有用
    # Refine knot vectors
    operations.refine_knotvector(surf, [3, 0])

    # Set visualization component
    vis_comp = VisMPL.VisSurface()
    surf.vis = vis_comp

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

instru = 'ME'
box_dex = 0
data_x,data_y,data_z,data_zerr = read_data(instru,box_dex)
a_bound,b_bound = get_bound(instru)

#-------tem for wrong zerr-------
data_zerr = np.sqrt(np.abs(data_z))
data_z = np.abs(data_z)

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
# data_z = data_z - factor
data_zerr = np.ones_like(data_zerr/np.max(data_zerr))
print(data_z)
print(np.max(data_z),np.min(data_z),np.mean(data_z))
print(data_zerr)

p_size_u = 200
p_size_v = 200
ctrl_x, ctrl_y = get_rectangle_divd(a_bound,b_bound,p_size_u,p_size_v)
p_ctrlpts, p_weights = get_ctrlpts(data_x,data_y,data_z,data_zerr,ctrl_x,ctrl_y)
print(np.max(p_ctrlpts[:,2]),np.min(p_ctrlpts[:,2]),np.mean(p_ctrlpts[:,2]))
p_degree_u = 20
p_degree_v = 20
psf_render(instru,box_dex,p_ctrlpts,p_weights,p_size_u, p_size_v, p_degree_u, p_degree_v)

psf = np.load('%s_psf_box%s.npy' % (instru, str(box_dex)))
print(np.max(psf[:,2]),np.min(psf[:,2]),np.mean(psf[:,2]))
print(psf)

