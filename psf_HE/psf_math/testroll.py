import numpy as np
import Quat as quat
from math import *
import time
def quat2mtx(qu):
    qu =list(qu)
    qu = qu+[sqrt(fabs(1-qu[0]*qu[0]-qu[1]*qu[1]-qu[2]*qu[2]))]
    q = quat.Quat(quat.normalize(qu))
    dcm_j_b = q.transform.T
    return dcm_j_b
def quat2delta2(mtx,dcm_b_f,xyz):
    #dcm_j_b =mtx 
    #vec_f = np.dot(dcm_b_f ,np.dot(dcm_j_b,xyz))
    vec_f = [mtx[1,0]*xyz[0]+mtx[1,1]*xyz[1]+mtx[1,2]*xyz[2],mtx[2,0]*xyz[0]+mtx[2,1]*xyz[1]+mtx[2,2]*xyz[2],mtx[0,0]*xyz[0]+mtx[0,1]*xyz[1]+mtx[0,2]*xyz[2]]
    alpha = atan(vec_f[0]/vec_f[2])*180/pi
    beta = atan(vec_f[1]/vec_f[2])*180/pi
    return alpha,beta

def quat2box(dcm_bx_to_b1,dcm_b_f, mtx, xyz):
    dcm_j_b =mtx
    vec_f = np.dot(dcm_bx_to_b1,np.dot(dcm_b_f,np.dot(dcm_j_b,xyz)))
    alpha = atan(vec_f[0]/vec_f[2])*180/pi
    beta = atan(vec_f[1]/vec_f[2])*180/pi
    return alpha,beta,vec_f[0],vec_f[1],vec_f[2]

def quat2delta2_rot(mtx,dcm_b_f,xyz,dcm_b_f_rot):
    dcm_j_b =mtx 
    vec_f = np.dot(dcm_b_f ,np.dot(dcm_b_f_rot,np.dot(dcm_j_b,xyz)))
    #print np.dot(dcm_b_f,dcm_b_f_rot)
    #vec_f = np.dot(dcm_b_f ,np.dot(dcm_j_b,xyz))
    #vec_f = [mtx[1,0]*xyz[0]+mtx[1,1]*xyz[1]+mtx[1,2]*xyz[2],mtx[2,0]*xyz[0]+mtx[2,1]*xyz[1]+mtx[2,2]*xyz[2],mtx[0,0]*xyz[0]+mtx[0,1]*xyz[1]+mtx[0,2]*xyz[2]]
    alpha = atan(vec_f[0]/vec_f[2])*180/pi
    beta = atan(vec_f[1]/vec_f[2])*180/pi
    return alpha,beta,vec_f[0],vec_f[1],vec_f[2]
def quat2delta(qu,ra,dec,roll=0):
#b = [-2.1583884e-01,4.9811391e-01,-4.0725711e-01,0.7345]
    qu =list(qu)
    #qu = qu+[sqrt(1-np.dot(qu,qu))]
    qu = qu+[sqrt(fabs(1-qu[0]*qu[0]-qu[1]*qu[1]-qu[2]*qu[2]))]
    #b=quat.normalize(qu)
    q = quat.Quat(quat.normalize(qu))
    #theta0 = 0
    angle = quat.Quat([0,0,0])
    [x,y,z] = sph2cart(ra*pi/180,dec*pi/180,1)
    dcm_j_b = q.transform.T
    #if roll==0: 
    #dcm1 = np.array([[0,1,0],[0,0,1],[1,0,0]])
    #if roll==60: 
    #    dcm1 = np.array([0,0.5,0.866,0,-0.866,0.5,1,0,0]).reshape(3,3)
    #if roll==-60: 
    #    dcm1 = np.array([0,0.5,-0.866,0,0.866,0.5,1,0,0]).reshape(3,3)
    #else:
    #    dcm1 = np.array([0,1,0,0,0,1,1,0,0]).reshape(3,3)
    #dcm2 = angle.transform.T
    #dcm_b_f = np.dot(dcm2,dcm1)
    dcm_b_f = np.array([[0,1,0],[0,0,1],[1,0,0]])
    vec_f = np.dot(dcm_b_f ,np.dot(dcm_j_b,[x,y,z]))
    alpha = atan(vec_f[0]/vec_f[2])*180/pi
    beta = atan(vec_f[1]/vec_f[2])*180/pi
    return alpha,beta
def sph2cart(az, el, r):
    rcos_theta = r * cos(el)
    x = rcos_theta * cos(az)
    y = rcos_theta * sin(az)
    z = r * sin(el)
    return x, y, z

def cart2sph(x, y, z):
    r=np.sqrt(x**2+y**2+z**2)
    lat = np.arcsin(z/r)*180/np.pi
    lon = np.arctan2(y,x)*180/np.pi
    tplon = lon if lon>0 else lon+360
    return lat,tplon

if __name__=="__main__":
    print("roll test")
    quat_list = np.loadtxt("/hxmt/work/jjin/forpanyy_forjguan_forychen/Output/G01HE0061ckli002/HOSIM_SatAttitude_G01HE0061ckli002_v0_00_new.dat")
    ra =  270
    dec = -30
    qu2 = [-2.1583884e-01,   4.9811391e-01,  -4.0725711e-01]
    a,b = quat2delta(qu2,-0.5666,-0.265435,0)
    print(qu2, a,b)
    start =time.clock()
    for i in range(50):
        for s_step in range(0,7639):
            qu3 = quat_list[s_step,:]
            delta_ra0,delta_dec0 = quat2delta(qu3,ra,dec,0)
        #    print delta_ra0,delta_dec0
    end = time.clock()
    print("time cost: ",end-start)

