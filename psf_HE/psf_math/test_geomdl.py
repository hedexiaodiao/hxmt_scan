#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: test_geomdl.py
# @Author: luoqi
# @Institution: Institute of High Energy Physics, CAS
# @E-mail: woshiluoqi123@outlook.com, luoqi@ihep.ac.cn
# @Site: 
# @Time: 12月 01, 2021
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

###需要有p_size_u >= p_degree_u+1, p_size_v >= p_degree_v+1, 一个是插值点, 一个是基元幂次数
# Dimensions of the control points grid
p_size_u = 4
p_size_v = 3

# Control points in u-row order
###控制点x、y坐标走向尽量按照等距先x变动再y变动的顺序
p_ctrlpts = [[0, 0, 0], [1, 0, 3], [2, 0, 0], [3, 0, 0],
             [0, 1, 0.2], [1, 1, 0.5], [2, 1, 0], [3, 1, -1],
             [0, 2, -0.6], [1, 2, 0], [2, 2, 3], [3, 2, 0]]

# Weights vector
p_weights = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

# Degrees
p_degree_u = 3
p_degree_v = 2


#
# Prepare data for import
#

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
surf.render()
