import numpy as np
import sys
import testroll###导入的旋转用模块
import Quat as quat ###python标准姿态四元数处理文件

def quat2mtx(qu):
    qu =list(qu)
    qu = qu+[np.sqrt(np.abs(1-qu[0]*qu[0]-qu[1]*qu[1]-qu[2]*qu[2]))]
    q = quat.Quat(quat.normalize(qu))
    dcm_j_b = q.transform.T
    return dcm_j_b

def quat2delta2_rot(dcm_b_f,dcm_j_b,xyz):
    vec_f = np.dot(dcm_b_f ,np.dot(dcm_j_b,xyz))
    alpha = np.atan(vec_f[0]/vec_f[2])*180/np.pi
    beta = np.atan(vec_f[1]/vec_f[2])*180/np.pi
    return alpha,beta,vec_f[0],vec_f[1],vec_f[2]

instr = "HE"

if instr == "HE":
    Talfa_bound=5.7 #means long_bound
    Tbeta_bound=1.1 #means short_bound
    dcm_b_f = np.array([[0,1,0],[0,0,1],[1,0,0]])
elif instr == "ME":
    Talfa_bound=4.0 #means long_bound
    Tbeta_bound=1.0 #means short_bound
    dcm_b_f = np.array([[0,0,-1],[0,1,0],[1,0,0]])
elif instr == "LE":
    Talfa_bound=6.0 #means long_bound
    Tbeta_bound=1.6 #means short_bound
    dcm_b_f = np.array([[0,0,1],[0,1,0],[-1,0,0]])
else:
    print("The instrment in config file is wrong.")
    sys.exit(1)

ra  = 262.991
dec = -33.834
quat1 = []###姿态四元数，q1,q2,q3
dcm_J2002payload = testroll.quat2mtx(quat1) ###从J2000到载荷
box_roll = [60,0,-60]###各机箱相对中间机箱旋转角度，采用了平面旋转
x, y, z = testroll.sph2cart(ra * np.pi / 180, dec * np.pi / 180, 1)
xyz = [x, y, z]
alpha, beta, x0, y0, z0 = quat2delta2_rot(dcm_b_f,dcm_J2002payload, xyz)
