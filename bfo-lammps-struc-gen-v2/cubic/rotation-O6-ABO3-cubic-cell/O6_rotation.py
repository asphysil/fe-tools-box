
import numpy as np
from gen_rotations import *

def O6_rotation_ABO3(nuo, nufe, axis, alpha, loc_axis_ordered, loc_noatm, atms_bo, o_atmu):

    o_atms_rot = np.zeros((nuo,3))

    for i in range(0, nufe):
        a1 = loc_axis_ordered[i, 1,:]
        a2 = loc_axis_ordered[i, 2,:]
        a3 = loc_axis_ordered[i, 3,:]
    
        trans_vec = loc_axis_ordered[i, 0,:]
       # print('trans',trans_vec)
    
        loc_no = loc_noatm[i]
    
        #vec_111 = a1 + a2 + a3
        # rotation w.r.t x-axis
        #vec_100=a1
        #vec_norm_100 = np.linalg.norm(vec_100)
        #u_vec100 = vec_100/vec_norm_100
        # rotation w.r.t y-axis
        #vec_010=a2 
        #--------------------Rotation matrix--------------
        if axis==1:
          rot_inphase = rotation_about_x(alpha)
          rot_antiphase = rotation_about_x(-alpha)
        elif axis==2:
            rot_inphase = rotation_about_y(alpha)
            rot_antiphase = rotation_about_y(-alpha)
        else:
            rot_inphase = rotation_about_z(alpha)
            rot_antiphase = rotation_about_z(-alpha)
    
    #print(fe_atmu)
        #print(atms_bo)
    
        if int(atms_bo[i][6]) == 1:
            fe=atms_bo[i,7:10]
            for j in range(0, loc_no):
                mm = int(atms_bo[i][j])
                vo =o_atmu[mm]
                vo_trans = vo -fe # trans_vec
    
                rvo = rotation_vec(rot_inphase, vo_trans)
                rvo = rvo + fe #trans_vec
                o_atms_rot[mm] =rvo
        else:
            fe=atms_bo[i,7:10]
            for j in range(0, loc_no):
                mm = int(atms_bo[i][j])
                vo =o_atmu[mm]
                vo_trans = vo - fe #trans_vec
    
                rvo = rotation_vec(rot_antiphase, vo_trans)
                rvo = rvo + fe #trans_vec
                o_atms_rot[mm] = rvo        
    return o_atms_rot
