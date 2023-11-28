
import numpy as np 


#q=dict()
#mass=dict()

def lammps_bv(natms, latt_sup, a_atms, b_atms, o_atms):

    nsa = natms[0]
    nsb = natms[1]
    nso = natms[2]
    # lammps data files
    flm = open("data.PMWO", 'w')
    flm.write('# LAMMPS ABO3 data file' + "\n\n")
    flm.write('  ' + '{0:d}'.format(nsa+nsb+nso) + ' ' + 'atoms' + "\n")
    flm.write('4' + '  ' + 'atom' + ' ' + 'types' + "\n\n")
    flm.write('0.0' + '    ' + str(latt_sup[0][0]) + '  ' + 'xlo' + ' ' + 'xhi' + "\n")
    flm.write('0.0' + '    ' + str(latt_sup[1][1]) + '  ' + 'ylo' + ' ' + 'yhi' + "\n")
    flm.write('0.0' + '    ' + str(latt_sup[2][2]) + '  ' + 'zlo' + ' ' + 'zhi' + "\n\n")
    flm.write('Masses' + "\n\n")
    a_mass = 207.2   # Pb
    b1_mass = 24.305  # Mg
    b2_mass = 183.84  # W
    o_mass = 15.99940 # O 
    flm.write('1' + '   ' +  str(a_mass) + "\n")
    flm.write('2' + '   ' +  str(b1_mass) + "\n")
    flm.write('3' + '   ' +  str(b2_mass) + "\n")
    flm.write('4' + '   ' +  str(o_mass) + "\n\n\n")
    flm.write('Atoms' +"\n\n")
    
    dl_bi = 0.00
    dl_fe = 0.00
    #
    ncount = 0
    nangle = 1
    
    a_chg = 1.77706    # Pb
    b1_chg = 1.93753   # Mg
    b2_chg = 1.93753   # W
    o_chg = -1.23820   # O
    
    ntype = 1
    for i in range(0, nsa):
        x = a_atms[i][0] 
        y = a_atms[i][1] 
        z = a_atms[i][2] 
        xqe = x 
        yqe = y 
        zqe = z 
        #fout_qe2.write('Bi   ' + '{0:12.8f} {1:12.8f} {2:12.8f}'.format(xqe/lx, yqe/ly, zqe/lz) + '  0  0  0' + "\n")
        ncount = ncount + 1
        flm.write('{0:d} {1:d} {2:d} {3:12.8f} {4:12.8f} {5:12.8f} {6:12.8f}'.format(ncount, nangle, ntype, a_chg, xqe, yqe, zqe) + "\n")
       
    for i in range(0, nsb):
        m = int(b_atms[i][0])
        x = b_atms[i][1]
        y = b_atms[i][2]
        z = b_atms[i][3]
        xqe = x 
        yqe = y 
        zqe = z 
        ncount = ncount + 1
        if m==1:
            ntype = 2
            flm.write('{0:d} {1:d} {2:d} {3:12.8f} {4:12.8f} {5:12.8f} {6:12.8f}'.format(ncount, nangle, ntype, b1_chg, xqe, yqe, zqe) + "\n")
        else:
            ntype = 3
            flm.write('{0:d} {1:d} {2:d} {3:12.8f} {4:12.8f} {5:12.8f} {6:12.8f}'.format(ncount, nangle, ntype, b2_chg, xqe, yqe, zqe) + "\n")
       
       
    ntype = 4    
    for i in range(0, nso):
        x = o_atms[i][0]
        y = o_atms[i][1]
        z = o_atms[i][2]
        ncount = ncount + 1
        flm.write('{0:d} {1:d} {2:d} {3:12.8f} {4:12.8f} {5:12.8f} {6:12.8f}'.format(ncount, nangle, ntype, o_chg, x, y, z) + "\n")
    