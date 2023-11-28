import numpy as np 

def sd_rombohedral_oxygen(n,  atms_in):
    
    d=0.01
    atms_out = np.zeros((n, 4))

    for i in range(0, n):
        x = atms_in[i][1]
        y = atms_in[i][2]
        z = atms_in[i][3] 

        xqe = x  + ((np.random.rand()-1.0)*d)
        yqe = y + ((np.random.rand()-1.0)*d)
        zqe = z + ((np.random.rand()-1.0)*d)

        atms_out[i][0] = atms_in[i][0]
        atms_out[i][1] = xqe 
        atms_out[i][2] = yqe 
        atms_out[i][3] = zqe
    return atms_out


def sd_rombohedral(n,  d111, latt, atms_in):

    atms_out = np.zeros((n, 3))

    # v1 = latt[0,:]
    # v2 = latt[1,:]
    # v3 = latt[2,:]

    # e1 = v1/np.linalg.norm(v1)
    # e2 = v2/np.linalg.norm(v2)
    # e3 = v3/np.linalg.norm(v3)
    # #print(e2)

    # r111 = e1+e2+e3

    # er = r111/np.linalg.norm(r111) 
    
    # vec_disp111 =  er*d111

    # e1_d = np.dot(vec_disp111, e1)*e1
    # e2_d = np.dot(vec_disp111, e2)*e2
    # e3_d = np.dot(vec_disp111, e3)*e3

    # #print(e1_d, e2_d, e3_d)

    # disp_vec = e1_d + e2_d + e3_d 

    # #print(np.linalg.norm(disp_vec), disp_vec, np.dot(disp_vec, e3))

    for i in range(0, n):
        x = atms_in[i][0]
        y = atms_in[i][1]
        z = atms_in[i][2] 

        xqe = x  + d111/np.sqrt(3) #disp_vec[0]
        yqe = y  + d111/np.sqrt(3) # disp_vec[1]
        zqe = z +  d111/np.sqrt(3) # disp_vec[2]

        atms_out[i][0] = xqe 
        atms_out[i][1] = yqe 
        atms_out[i][2] = zqe
    return atms_out

def sd_tetragonal(n,  d001, atms_in):

    atms_out = np.zeros((n, 4))

    for i in range(0, n):
        x = atms_in[i][1]
        y = atms_in[i][2]
        z = atms_in[i][3] 

        xqe = x 
        yqe = y
        zqe = z + d001

        atms_out[i][0] = atms_in[i][0]
        atms_out[i][1] = xqe 
        atms_out[i][2] = yqe 
        atms_out[i][3] = zqe
    return atms_out