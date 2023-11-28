import numpy as np 


def BO_1stdist_pbc(nn, nb, no, latt_vec, dbo, atmb, atmo):
    m1 = 0
    m2 = 0
    min_dist = np.zeros(no*27, dtype=np.float)
    o_index = np.zeros(no*27, dtype=np.int)
   
    atms_bo = np.ones((nb, 7))*(-1)
    vec_n = np.zeros(3)


    for i in range(0, nb):
        atms_bo[i][6] = atmb[i][0]
        
        vb = atmb[i,1:4]
        mm = 0 
        for j in range(0, no):
            vo = atmo[j]
            for ix in range(-1, 2):
                for iy in range(-1, 2):
                    for iz in range(-1, 2):
                        vec_n[0] = ix 
                        vec_n[1] = iy 
                        vec_n[2] = iz 
                        T = np.matmul(np.transpose(latt_vec), vec_n)
                    #print(T)
                        w = vo + T  
                        diff = vb-w 
                        min_dist[mm] = np.linalg.norm(diff)
                        o_index[mm] = j
                        mm = mm + 1

        for k in range(0, nn):
            m = np.argmin(min_dist)
            d = min_dist[m]
            if d <= dbo:
                atms_bo[i][k] = o_index[m]
                min_dist[m] = 10000.0
            else:
                pass
    
   # print(temp_atmbo, np.argwhere(temp_atmbo[0] == -1))

    #atms_bo = np.delete(temp_atmbo, np.argwhere(temp_atmbo[0] == -1))
    #print(atms_bo)

    return (atms_bo)

def BO_1stdist_nonpbc(nn, nb, no, dbo, atmb, atmo):
    m1 = 0
    m2 = 0
    min_dist = np.zeros(no, dtype=float)
   
    atms_bo = np.ones((nb, 10))*(-1)


    for i in range(0, nb):
        atms_bo[i][6] = atmb[i][0] # 
        atms_bo[i,7:10] = atmb[i][1:4]  # Fe-atomic coordinates
        vb = atmb[i,1:4]
        mm = 0 
        for j in range(0, no):
            vo = atmo[j]
            diff = vb-vo 
            min_dist[j] = np.linalg.norm(diff)

        for k in range(0, nn):
            m = np.argmin(min_dist)
            d = min_dist[m]
            if d <= dbo:
                atms_bo[i][k] = m
                min_dist[m] = 10000.0
            else:
                pass
    
   # print(temp_atmbo, np.argwhere(temp_atmbo[0] == -1))

    #atms_bo = np.delete(temp_atmbo, np.argwhere(temp_atmbo[0] == -1))
    #print(atms_bo)

    return (atms_bo)

def AB_1stdist_pbc(nn, na, nb, dab, latt_vec, atma, atmb):
    m1 = 0
    m2 = 0
    min_dist = np.zeros(na*27, dtype=float)
    a_index = np.zeros(na*27, dtype=float)
    a_coord = np.zeros((na*27,3), dtype=float)

#-------------------------------------------------
    coord_ab = np.zeros((nb,8,3), dtype=float)
    atms_ab = np.ones((nb, 8))*(-1)
    vec_n = np.zeros(3)


    for i in range(0, nb):
        vb = atmb[i,1:4]
        mm = 0 
   
        for j in range(0, na):
            va = atma[j]
            for ix in range(-1, 2):
                for iy in range(-1, 2):
                    for iz in range(-1, 2):
                        vec_n[0] = ix 
                        vec_n[1] = iy 
                        vec_n[2] = iz 
                        T = np.matmul(np.transpose(latt_vec), vec_n)
                    #print(T)
                        w = va + T  
                        diff = vb-w 
                        min_dist[mm] = np.linalg.norm(diff)
                        a_index[mm] = j
                        a_coord[mm] = w
                        mm = mm + 1

        for k in range(0, nn):
            m = np.argmin(min_dist)
            d = min_dist[m]
            if d <= dab:
                atms_ab[i][k] = a_index[m]
                coord_ab[i][k] = a_coord[m]
                min_dist[m] = 10000.0
            else:
                pass
    
   # print(temp_atmbo, np.argwhere(temp_atmbo[0] == -1))

    #atms_bo = np.delete(temp_atmbo, np.argwhere(temp_atmbo[0] == -1))
    #print(atms_bo)

    return (atms_ab, coord_ab)


def center_of_mass_O6(latt_vec, nn, nb, no, dbo, atmb, atmo):
    
    min_dist = np.zeros(27*no)
    coord_o = np.zeros((27*no, 3))
    vec_n = np.zeros(3)
    cmass = np.zeros((nb,3))
    
    b_atms_disp_O6 = np.zeros((nb,3))
    
    for i in range(0, nb):
        vb = atmb[i]
        vf = np.zeros(3)

        nn = 0 
        for j in range(0, no):
            vo = atmo[j]
            for ix in range(-1, 2):
                for iy in range(-1, 2):
                    for iz in range(-1, 2):
                        vec_n[0] = ix 
                        vec_n[1] = iy 
                        vec_n[2] = iz 
                        T = np.matmul(np.transpose(latt_vec), vec_n)
                    #print(T)
                        w = vo + T  
                        diff = vb-w 
                        min_dist[nn] = np.linalg.norm(diff)
                        coord_o[nn] = w 
                        nn = nn + 1
                    
            #print(min_dist)
        for k in range(0,6):
            mm = np.argmin(min_dist)
                #print(coord_o[mm])
            vf = coord_o[mm] + vf 
            min_dist[mm] = 1000
  
        cmass[i]= vf/6.0 


    for i in range(0, nb):
        b_atms_disp_O6[i] = atmb[i]-cmass[i]

    #np.set_printoptions(precision=3)
    #print(cmass)
    print(np.around(b_atms_disp_O6, decimals=3))
    return 

def center_of_mass_O12(latt_vec, nn, nb, no, dbo, atma, atmo):
    
    min_dist = np.zeros(27*no)
    coord_o = np.zeros((27*no, 3))
    vec_n = np.zeros(3)
    cmass = np.zeros((nb,3))
    
    a_atms_disp_O12 = np.zeros((nb,3))
    
    for i in range(0, nb):
        vb = atma[i]
        vf = np.zeros(3)

        nn = 0 
        for j in range(0, no):
            vo = atmo[j]
            for ix in range(-1, 2):
                for iy in range(-1, 2):
                    for iz in range(-1, 2):
                        vec_n[0] = ix 
                        vec_n[1] = iy 
                        vec_n[2] = iz 
                        T = np.matmul(np.transpose(latt_vec), vec_n)
                    #print(T)
                        w = vo + T  
                        diff = vb-w 
                        min_dist[nn] = np.linalg.norm(diff)
                        coord_o[nn] = w 
                        nn = nn + 1
                    
            #print(min_dist)
        for k in range(0,12):
            mm = np.argmin(min_dist)
                #print(coord_o[mm])
            vf = coord_o[mm] + vf 
            min_dist[mm] = 1000
  
        cmass[i]= vf/12.0 


    for i in range(0, nb):
        a_atms_disp_O12[i] = atma[i]-cmass[i]

    #np.set_printoptions(precision=3)
    print(cmass)
    print(np.around(a_atms_disp_O12, decimals=5))
    return 

def cal_local_axis(nn, nb, atms_coord):

    local_axis = np.zeros((nb, 4, 3))
    min_dist = np.zeros(7)

    for i in range(0, nb):
        local_a_atms = np.zeros((8,3))
        local_a_atms = atms_coord[i]

        #my_array[:, [2, 0]] = my_array[:, [0, 2]]

        print(local_a_atms)
        local_a_atms = local_a_atms[local_a_atms[:, 1].argsort()]
        local_a_atms = local_a_atms[local_a_atms[:, 2].argsort()]
        local_a_atms = local_a_atms[local_a_atms[:, 0].argsort()]

        local_axis[i,0:1]=local_a_atms[0]

        v1 = local_a_atms[0] 
        print('==========')
        for k1 in range(1, nn):
            v2 = local_a_atms[k1]
            print(v2) 
            diff = v1-v2
            min_dist[k1-1] = np.linalg.norm(diff)
        #print(local_a_atms)
        for k2 in range(0,3):
            mm = np.argmin(min_dist)
            local_axis[i, k2+1] = local_a_atms[mm+1]
            min_dist[mm] = 10000.0

            


        #print(local_a_atms)


        #print(local_a_atms)
    print(local_axis) 

    return local_axis



def local_axis_arrange(nn_bo, nb, local_axis, atms_bo):
    

    loc_axis_ordered = np.zeros((nb,4, 3))

    loc_axis_natms = np.zeros(nb, dtype=int)

    for i in range(0, nb):
        a1=local_axis[i, 3:4]-local_axis[i,0:1]
        a2=local_axis[i, 1:2]-local_axis[i,0:1]
        a3=local_axis[i, 2:3]-local_axis[i,0:1]

        loc_axis_ordered[i,0,0:3] = local_axis[i,0:1] 
        loc_axis_ordered[i,1,0:3] = a1
        loc_axis_ordered[i,2,0:3] = a2
        loc_axis_ordered[i,3,0:3] = a3
    
        #print(local_latt)
        #print(atms_bo)
        nn = 0
        for j in range(0, nn_bo):
            mm = int(atms_bo[i][j])
            if mm != -1:
               nn = nn + 1
               #print(o_frac)
        loc_axis_natms[i]=nn


    #print(loc_axis_and_oatms_coord)
    #print(loc_axis_natms)
    return (loc_axis_natms, loc_axis_ordered)


def local_axis_and_oatms(nn_bo, nb, local_axis, atms_bo, atms_o):

    local_latt = np.zeros((3,3))

    loc_axis_and_oatms_coord = np.zeros((nb, 9, 3))

    loc_axis_natms = np.zeros(nb, dtype=np.int)

    for i in range(0, nb):
        a1=local_axis[i, 3:4]-local_axis[i,0:1]
        a2=local_axis[i, 1:2]-local_axis[i,0:1]
        a3=local_axis[i, 2:3]-local_axis[i,0:1]

        loc_axis_and_oatms_coord[i,0,0:3] = a1
        loc_axis_and_oatms_coord[i,1,0:3] = a2
        loc_axis_and_oatms_coord[i,2,0:3] = a3
    

        local_latt[0:,] = a1
        local_latt[1:,] = a2
        local_latt[2:,] = a3
        
        #print(local_latt)
        #print(atms_bo)
        nn = 0
        for j in range(0, nn_bo):
            mm = np.int(atms_bo[i][j])
            if mm != -1:
               vo =atms_o[mm] 
               vo_trans = vo-local_axis[i,0:1]

               o_frac = np.matmul(vo_trans, np.linalg.inv(local_latt))
               loc_axis_and_oatms_coord[i,3+nn,0:3] = o_frac
 
               nn = nn + 1
               #print(o_frac)
        loc_axis_natms[i]=nn


    #print(loc_axis_and_oatms_coord)
    #print(loc_axis_natms)
    return (loc_axis_natms, loc_axis_and_oatms_coord)