import numpy as np
 

fin_lammps='data-jiahao.BFO'

fout_lammps='bfo-02.data'


f1=open(fout_lammps, 'r')
n=len(f1.readlines())
print(n)
f1.close()

f1=open(fin_lammps, 'r')
n_in=len(f1.readlines())
print(n_in)
f1.close()


l1 = np.genfromtxt(fout_lammps, skip_header=1, skip_footer=n-14, dtype=str)
ntypes=np.genfromtxt(fout_lammps, skip_header=3, skip_footer=n-15, dtype=str)
print(ntypes)
nangles=np.genfromtxt(fout_lammps, skip_header=4, skip_footer=n-16, dtype=str)
print(nangles)
angle_types=np.genfromtxt(fout_lammps, skip_header=5, skip_footer=n-17, dtype=str)
print(angle_types)
#
a1 = np.genfromtxt(fout_lammps, skip_header=7, skip_footer=n-18, dtype=str)
#
a2 = np.genfromtxt(fout_lammps, skip_header=8, skip_footer=n-19, dtype=str)
#
a3 = np.genfromtxt(fout_lammps, skip_header=9, skip_footer=n-20, dtype=str)
#
#l3 = np.genfromtxt(fout_lammps, skip_header=1, skip_footer=n-16, dtype=str)
#print(l3)
#print(l1)
#print(a1, a2, a3)
latt_sup=np.zeros((3,3))

latt_sup[0,0] = float(a1[1]) - float(a1[0])  
latt_sup[1,1] = float(a2[1]) - float(a2[0])
latt_sup[2,2] = float(a3[1]) - float(a3[0])

natms = int(l1[0])
nangles = int(nangles[0])
ntypes = int(ntypes[0])
angle_types = int(angle_types[0])

print(ntypes)
print(natms)
print(nangles)

# output data
m = n-natms-26
data = np.genfromtxt(fout_lammps, usecols=(0,1,3,4,5,6,7), skip_header=22, skip_footer=m)
print(data)

# Angle data
p=n-nangles
angle_data=np.genfromtxt(fout_lammps, skip_header=p) #, skip_footer=m)
#print(angle_data)

# Arrange data
data[:,[0,6]] = data[:,[6,0]]
atms = data[data[:,6].argsort()]
atms[:,[0,6]] = atms[:,[6,0]]

# orginal data 
m = n_in-natms-22
atms_org = np.genfromtxt(fin_lammps, skip_header=20, skip_footer=m)
#print(atms_org)

def write_lammps_bv_input(natms_types, nangles_types, natms, nangles, latt_sup, atms_org, atms, angle):

    # lammps data files
    flm = open("data-out.BFO", 'w')

    flm.write('# LAMMPS ABO3 data file' + "\n\n")
    flm.write('  ' + '{0:d}'.format(natms) + ' ' + 'atoms' + "\n")
    flm.write('  ' + '{0:d}'.format(nangles) + ' ' + 'angles' + "\n")
    flm.write(str(natms_types) + '  ' + 'atom' + ' ' + 'types' + "\n\n")
    flm.write(str(nangles_types) + '  ' + 'angle' + ' ' + 'types' + "\n\n")
    flm.write('0.0' + '    ' + str(latt_sup[0][0]) + '  ' + 'xlo' + ' ' + 'xhi' + "\n")
    flm.write('0.0' + '    ' + str(latt_sup[1][1]) + '  ' + 'ylo' + ' ' + 'yhi' + "\n")
    flm.write('0.0' + '    ' + str(latt_sup[2][2]) + '  ' + 'zlo' + ' ' + 'zhi' + "\n\n")
    flm.write('Masses' + "\n\n")
    a_mass = 208.9804   # Bi
    b_mass = 55.845  # Fe
    o_mass = 15.99940 # O 
    flm.write('1' + '   ' +  str(a_mass) + "\n")
    flm.write('2' + '   ' +  str(b_mass) + "\n")
    flm.write('3' + '   ' +  str(o_mass) + "\n\n")
    flm.write('Atoms' +"\n\n")

    #
    ncount = 0
    nangle = 1

    for i in range(0, natms):
        flm.write(' {:d} '.format(int(atms[i,0])))
        flm.write(' {:d} '.format(int(atms[i,1])))
        flm.write(' {:d} '.format(int(atms_org[i,2])))
        for j in range(2,7):
            flm.write(' {:12.8f} '.format(atms[i,j]))

        for j in range(0,3):
            flm.write(' {:d} '.format(int(atms_org[i,j+8])))
        flm.write("\n")


    flm.write("\n")
    flm.write('Angles' +"\n\n")
    for i in range(0, nangles):
          for j in range(0, 5):
              flm.write(' {:d}  '.format(int(angle[i,j])))
          flm.write("\n")

    flm.close()

write_lammps_bv_input(ntypes, angle_types, natms, nangles, latt_sup, atms_org, atms, angle_data)

