import numpy as np 

def CtypeAFM_order_zaxis(nx, ny, nz, nb, latt, atms):
    nn=4
    nyz=ny*nz  # Number of atom in yz-plane
    nxy=nx*ny  # Number of atom in xy-plane
    nxz=nx*nz # Number of atoms in xz-plane 

    dmin = min_dist(nb, atms) + 1.0 # minimum distance between two atoms
    #print(dmin)

    dist = np.zeros(nxy)  # plane perpendicular to z-axis 

    atms_order = np.zeros((nb,5)) #
    dist_pbc = np.zeros(27)
    v2 = np.zeros(3)
    vec_n = np.zeros(3)

    for i in range(0, nb):
        atms_order[i][0] = i
        atms_order[i][1] = 3
        atms_order[i][2] = atms[i][0]
        atms_order[i][3] = atms[i][1]
        atms_order[i][4] = atms[i][2]

    #atms_order[:, [4, 4]] = atms_order[:, [4, 4]] # Column swape
    atms_order=atms_order[atms_order[:,4].argsort()]
    #print(atms_order)
    atms_order_temp=atms_order[0:nxy,:] # copy 
    #print(" ")
   # print(atms_order_temp)
    atms_order_temp[:,1]=1.0

    # G-type-AFM  order at the bottom layer 
    w=np.zeros(3)
    for i in range(0, nxy):
        v1 = atms_order_temp[i,2:5]
        
        if int(atms_order_temp[i][1])==1:
            print(i)
            
            for j in range(0, nxy):
                v2=atms_order_temp[j,2:5]
                m = 0
                for lx in range(-1, 2):
                    for ly in range(-1,2):
                        for lz in range(-1,2):
                            vec_n[0] = lx 
                            vec_n[1] = ly 
                            vec_n[2] = lz 
                            T = np.matmul(np.transpose(latt), vec_n)
                            #print(T) 
                            w[0] = v2[0] + T[0]
                            w[1] = v2[1] + T[1]
                            w[2] = v2[2] + T[2]
                            dist_pbc[m] = np.linalg.norm(v1-w)
                            m = m + 1
                dist[j] = np.min(dist_pbc) 
            m1 = np.argmin(dist)
            dist[m1] = 10000.0
            #print(dist)
            for k in range(0, nn):
                m1 = np.argmin(dist)
                d = dist[m1]
                if d <= dmin:
                   # print(int(atms_order_temp[m][1]))
                    dist[m1] = 10000.0
                    if int(atms_order_temp[m1][1]) == 1:
                         #print('ok')
                         atms_order_temp[m1][1] = 2.0
                         #print(m1)
        else:
            pass
    
    atms_order[0:nxy, :]=atms_order_temp[0:nxy,:] # back to orginal array
    #print(atms_order)
    #print("ok")

    # To get atoms along the x-axis and yz-plane
    atms_order[:, [4, 2]] = atms_order[:, [2, 4]] # Column swapped x and z
    #print(atms_order)
    #print("----")
    atms_order=atms_order[atms_order[:,4].argsort()]
    #print(atms_order)
    n1=0
    n2=nyz 
    for i in range(0, nx): # number of layer along x-axis
        atms_order_temp=atms_order[n1:n2,:] # Number of atoms in yz plane 
        atms_order_temp[:, [4, 3]] = atms_order_temp[:, [3, 4]] # Column swapped y and z
        atms_order_temp=atms_order_temp[atms_order_temp[:,4].argsort()]  
        m1=0
        m2=nz 
        for i in range(0, ny): # number of layer along y-axis
            atms=atms_order_temp[m1:m2,:] # Number of atoms along z-direction
            d=atms[0:nz, 1]  
            l=np.min(d)
            atms_order_temp[m1:m2,1]=l
            m1=m2
            m2 += nz 

        atms_order_temp[:, [4, 3]] = atms_order_temp[:, [3, 4]] # Column swapped y and z 
        atms_order[n1:n2,:]=atms_order_temp # put back to orginal array 

        n1=n2 
        n2 +=nyz 

        
    atms_order[:, [4, 2]] = atms_order[:, [2, 4]] # Column swapped x and z 
    
    atms_order[:, [4, 0]] = atms_order[:, [0, 4]] # Column swapped  0 and z 
    atms_order=atms_order[atms_order[:,4].argsort()]    # back original order
    atms_order[:, [4, 0]] = atms_order[:, [0, 4]] # Column swapped back 
    return atms_order



def CtypeAFM_order_xaxis(nx, ny, nz, nb, latt, atms):
    nn=4
    nyz=ny*nz 
    nxy=nx*ny 
    nxz=nx*nz

    dmin = min_dist(nb, atms) + 1.0 # minimum distance between two atoms
    #print(dmin)

    dist = np.zeros(nyz) # plane perpendicular to x-axis 

    atms_order = np.zeros((nb,5)) #
    dist_pbc = np.zeros(27)
    v2 = np.zeros(3)
    vec_n = np.zeros(3)

    for i in range(0, nb):
        atms_order[i][0] = i
        atms_order[i][1] = 3
        atms_order[i][2] = atms[i][0]
        atms_order[i][3] = atms[i][1]
        atms_order[i][4] = atms[i][2]

    atms_order[:, [4, 2]] = atms_order[:, [2, 4]] # Column swapped  z and x
    atms_order=atms_order[atms_order[:,4].argsort()]
    #print(atms_order)
    atms_order_temp=atms_order[0:nyz,:]
    #print(" ")
    #print(atms_order_temp)
    atms_order_temp[:,1]=1.0

    # G-type-AFM  order at the bottom layer 
    w=np.zeros(3)
    for i in range(0, nyz):
        v1 = atms_order_temp[i,2:5]
        
        if int(atms_order_temp[i][1])==1:
           # print(i)
            
            for j in range(0, nyz):
                v2=atms_order_temp[j,2:5]
                m = 0
                for lx in range(-1, 2):
                    for ly in range(-1,2):
                        for lz in range(-1,2):
                            vec_n[0] = lx 
                            vec_n[1] = ly 
                            vec_n[2] = lz 
                            T = np.matmul(np.transpose(latt), vec_n)
                            #print(T) 
                            w[0] = v2[0] + T[0]
                            w[1] = v2[1] + T[1]
                            w[2] = v2[2] + T[2]
                            dist_pbc[m] = np.linalg.norm(v1-w)
                            m = m + 1
                dist[j] = np.min(dist_pbc) 
            m1 = np.argmin(dist)
            dist[m1] = 10000.0
            #print(dist)
            for k in range(0, nn):
                m1 = np.argmin(dist)
                d = dist[m1]
                if d <= dmin:
                   # print(int(atms_order_temp[m][1]))
                    dist[m1] = 10000.0
                    if int(atms_order_temp[m1][1]) == 1:
                         #print('ok')
                         atms_order_temp[m1][1] = 2.0
                         #print(m1)
        else:
            pass
    
    atms_order[0:nyz, :]=atms_order_temp[0:nyz,:] # back to orginal array
    #print(atms_order)
    #print("ok")
    # back to orginal
    atms_order[:, [4, 2]] = atms_order[:, [2, 4]] # Column swapped x and z
    #print(atms_order)
    #print("----")

    #atms_order[:, [4, 3]] = atms_order[:, [3, 4]] # Column swapped y and z
    atms_order=atms_order[atms_order[:,4].argsort()]
    #print(atms_order)
    n1=0
    n2=nxy
    for i in range(0, nz): # number of layer along y-axis
        atms_order_temp=atms_order[n1:n2,:] # Number of atoms in xy plane 
        atms_order_temp[:, [4, 3]] = atms_order_temp[:, [3, 4]] # Column swapped x and z
        atms_order_temp=atms_order_temp[atms_order_temp[:,4].argsort()]  
        m1=0
        m2=nx  # Number of atom along x-axis
        for i in range(0, ny): # number of layer along y-axis
            atms=atms_order_temp[m1:m2,:] # Number of atoms along x-direction
            #print(atms)
            d=atms[0:nx, 1]  
            l=np.min(d)
            print(l)
            atms_order_temp[m1:m2,1]=l
            m1=m2
            m2 += nx

        atms_order_temp[:, [4, 3]] = atms_order_temp[:, [3, 4]] # Column swapped y and z 
        atms_order[n1:n2,:]=atms_order_temp # put back to orginal array 

        n1=n2 
        n2 +=nxy

        
    #atms_order[:, [4, 3]] = atms_order[:, [3, 4]] # Column swapped x and z 
    #print(atms_order)

    atms_order[:, [4, 0]] = atms_order[:, [0, 4]] # Column swapped  0 and z 
    atms_order=atms_order[atms_order[:,4].argsort()]    # back original order
    atms_order[:, [4, 0]] = atms_order[:, [0, 4]] # Column swapped back 
    return atms_order



# 
def CtypeAFM_order_yaxis(nx, ny, nz, nb, latt, atms):
    nn=4
    nyz=ny*nz 
    nxy=nx*ny 
    nxz=nx*nz

    dmin = min_dist(nb, atms) + 1.0 # minimum distance between two atoms
    #print(dmin)

    dist = np.zeros(nxz) # plane perpendicular to y-axis 

    atms_order = np.zeros((nb,5)) #
    dist_pbc = np.zeros(27)
    v2 = np.zeros(3)
    vec_n = np.zeros(3)

    for i in range(0, nb):
        atms_order[i][0] = i
        atms_order[i][1] = 3
        atms_order[i][2] = atms[i][0]
        atms_order[i][3] = atms[i][1]
        atms_order[i][4] = atms[i][2]

    atms_order[:, [4, 3]] = atms_order[:, [3, 4]] # Column swapped  z and y
    atms_order=atms_order[atms_order[:,4].argsort()]
    #print(atms_order)
    atms_order_temp=atms_order[0:nxz,:]
    #print(" ")
    #print(atms_order_temp)
    atms_order_temp[:,1]=1.0

    # G-type-AFM  order at the bottom layer 
    w=np.zeros(3)
    for i in range(0, nxz):
        v1 = atms_order_temp[i,2:5]
        
        if int(atms_order_temp[i][1])==1:
           # print(i)
            
            for j in range(0, nxz):
                v2=atms_order_temp[j,2:5]
                m = 0
                for lx in range(-1, 2):
                    for ly in range(-1,2):
                        for lz in range(-1,2):
                            vec_n[0] = lx 
                            vec_n[1] = ly 
                            vec_n[2] = lz 
                            T = np.matmul(np.transpose(latt), vec_n)
                            #print(T) 
                            w[0] = v2[0] + T[0]
                            w[1] = v2[1] + T[1]
                            w[2] = v2[2] + T[2]
                            dist_pbc[m] = np.linalg.norm(v1-w)
                            m = m + 1
                dist[j] = np.min(dist_pbc) 
            m1 = np.argmin(dist)
            dist[m1] = 10000.0
            #print(dist)
            for k in range(0, nn):
                m1 = np.argmin(dist)
                d = dist[m1]
                if d <= dmin:
                   # print(int(atms_order_temp[m][1]))
                    dist[m1] = 10000.0
                    if int(atms_order_temp[m1][1]) == 1:
                         #print('ok')
                         atms_order_temp[m1][1] = 2.0
                         #print(m1)
        else:
            pass
    
    atms_order[0:nxz, :]=atms_order_temp[0:nxz,:] # back to orginal array
    #print(atms_order)
    #print("ok")
    # back to orginal
    atms_order[:, [4, 3]] = atms_order[:, [3, 4]] # Column swapped y and z
    #print(atms_order)
    #print("----")

    #atms_order[:, [4, 3]] = atms_order[:, [3, 4]] # Column swapped y and z
    atms_order=atms_order[atms_order[:,4].argsort()]
    #print(atms_order)
    n1=0
    n2=nxy
    for i in range(0, nz): # number of layer along y-axis
        atms_order_temp=atms_order[n1:n2,:] # Number of atoms in xy plane 
        atms_order_temp[:, [4, 2]] = atms_order_temp[:, [2, 4]] # Column swapped x and z
        atms_order_temp=atms_order_temp[atms_order_temp[:,4].argsort()]  
        m1=0
        m2=ny 
        for i in range(0, nx): # number of layer along x-axis
            atms=atms_order_temp[m1:m2,:] # Number of atoms along z-direction
            #print(atms)
            d=atms[0:ny, 1]  
            l=np.min(d)
            print(l)
            atms_order_temp[m1:m2,1]=l
            m1=m2
            m2 += ny

        atms_order_temp[:, [4, 2]] = atms_order_temp[:, [2, 4]] # Column swapped y and z 
        atms_order[n1:n2,:]=atms_order_temp # put back to orginal array 

        n1=n2 
        n2 +=nxy

        
    #atms_order[:, [4, 3]] = atms_order[:, [3, 4]] # Column swapped x and z 
    #print(atms_order)

    atms_order[:, [4, 0]] = atms_order[:, [0, 4]] # Column swapped  0 and z 
    atms_order=atms_order[atms_order[:,4].argsort()]    # back original order
    atms_order[:, [4, 0]] = atms_order[:, [0, 4]] # Column swapped back 
    return atms_order


#-----------------------------------


def min_dist(n, atm):

    dist = np.zeros(n-1)
    
    for i in range(0,1):
        v1 = atm[i]
        m = 0
        for j in range(0, n):
            if j != i:
                #print(j)
                v2 = atm[j]
                dv = v1-v2
                dist[m] = np.linalg.norm(dv)
                m = m +1
                
            
    d = np.min(dist)
    return d