#!/usr/bin/python

import sys
from numpy import loadtxt, zeros,argsort, tanh
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

file_name='cell_P.dat'
x1, x2, x3, px, py, pz=loadtxt(file_name, unpack=True)

nx=8
ny=1
nz=1
nti=nx*ny*nz

if nti==len(x1):
    pass
else:
    print('dimension of the cell is not correct')
    sys.exit()

atm_info=zeros((len(x1),6))

for i in range(0, len(x1)):
    atm_info[i][0]=px[i]
    atm_info[i][1]=py[i]
    atm_info[i][2]=pz[i]
    atm_info[i][3]=x3[i]
    atm_info[i][4]=x2[i]
    atm_info[i][5]=x1[i]
    
atm_arrange_x=atm_info[atm_info[:,5].argsort()]

nyzatms=ny*nz

nxlayer=zeros(nx)
px=zeros(nx)
py=zeros(nx)
pz=zeros(nx)
c_ti=zeros(nx)

Ti=22

n=0
for i in range(0, nx):
    nxlayer[i]=i
    for j in range(0, nyzatms):
        px[i]=px[i]+atm_arrange_x[n][0]
        py[i]=py[i]+atm_arrange_x[n][1]
        pz[i]=pz[i]+atm_arrange_x[n][2]
        #print(atm_arrange_x[n][3], atm_arrange_x[n][6])
        c_ti[i]=c_ti[i]+ 1
        n=n+1

px=px/nyzatms
py=py/nyzatms
pz=pz/nyzatms
c_ti=c_ti/nyzatms

for i in range(0, nx):
    s='{:0.3}'.format(c_ti[i])
    c_ti[i]=s

def func_tanh(x, x0, zhi, pz0):
    return pz0*tanh((x-x0)/zhi)

plt.plot(nxlayer, px, 'o--', label='px')
plt.plot(nxlayer, py, 'H--', label='py')

plt.plot(nxlayer, pz, 'D--', label='pz')

pz_fit=zeros(nx)

pz_fit[0]=pz[1]
pz_fit[len(pz)-1]=pz[len(pz)-2]

pz_fit[1:len(pz)-1]=pz[1:len(pz)-1]

popt, pcov = curve_fit(func_tanh, nxlayer, pz_fit, p0=[4, 1.0, 0.68])

print(popt)

pz_fit=func_tanh(nxlayer, popt[0], popt[1], popt[2])

plt.plot(nxlayer, pz_fit, 'd--', label='pz_fit')

#plt.xticks(nxlayer, c_ti)
plt.tick_params(axis='y', which='both', labelleft=True, labelright=True)
plt.xlim([0,11])
plt.legend()
plt.ylabel('Polarization (C/m$^2$)')
plt.xlabel('Layer No')
plt.show()
