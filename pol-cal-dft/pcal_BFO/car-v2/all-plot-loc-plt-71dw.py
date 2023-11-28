import numpy as np
import time


import matplotlib.pyplot as plt



n = 1

nx=10
ny=2
nz=4
nuti=2

nti = (nx*ny*nz)*nuti

file_name='cell_P.dat'

#k = 501

for k in range(1, 2):
    mi = (k - 1) * nti
    mf = (n - k) * nti
   
 
    data = np.genfromtxt(file_name, skip_header = mi, skip_footer = mf)
 
    dt = 0.0005
    
    tn = 1000*dt*n
    
    atm_info = np.zeros((nti, 6))
    
    for l in range(0, 1):
        kk = l*nti
        for i in range(0, nti):
            m = kk + i
            atm_info[i][0] = data[m][3]
            atm_info[i][1] = data[m][4]
            atm_info[i][2] = data[m][5]
            atm_info[i][3] = data[m][1]
            atm_info[i][4] = data[m][2]
            atm_info[i][5] = data[m][0]
    
    
        ploc = atm_info[atm_info[:,5].argsort()]
        x = np.zeros(nti)
        y0 = np. zeros(nti)
        y1 = np. zeros(nti)
        y2 = np. zeros(nti)
        for i in range(0, nti):
            x[i]  =ploc[i][5]
            y0[i] =ploc[i][0]
            y1[i] =ploc[i][1]
            y2[i] =ploc[i][2]
    
    
    plt.plot(x,y0, 'o--', label = 'p-->[1-10]')

    plt.plot(x,y1, 'H--', label = 'p-->[110]')
    plt.plot(x,y2, 'D--', label = 'p-->[001]')
    #plt.plot(tn, p, label = 'p = $\sqrt(px^2 + py^2 + pz^2)')
    
    plt.ylabel('Polarization per unit cell (C/m$^{2}$)')
    plt.xlabel('Distance --> [001]')
    plt.legend(loc='best')
#    f='fig-'+ str(k) + '.png'
    #print(f)
#    plt.savefig(f)
    #plt.pause(0.5)
    #plt.close()
plt.show()
    #time.sleep(1.0)
