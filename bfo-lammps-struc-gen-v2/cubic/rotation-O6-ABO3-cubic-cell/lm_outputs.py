import numpy as np 


def siesta_output(nspecies, natms, atomic_num, atms_name, latt_vec, bi_atms, fe_atms,  o_atms):
    ntotal_atms = np.sum(natms)
    ntype_order = np.zeros(nspecies, dtype=np.int)

    a = latt_vec[0][0] # [1-10]
    b = latt_vec[1][1] # [110]
    c = latt_vec[2][2] # [001]

    fout_siesta=open('STRUC_16L.fdf', 'w')
    fout_siesta_spin=open('spin_input.siesta', 'w')

    fout_siesta.write('# structural and geometry variables' +"\n")
    fout_siesta.write('{0:s} {1:d}'.format('NumberOfAtoms', ntotal_atms) + "\n")
    fout_siesta.write('{0:s} {1:d}'.format('NumberOfSpecies', nspecies) + "\n")
    fout_siesta.write('{0:s}'.format('LatticeConstant       1.0 Ang')+"\n")
    fout_siesta.write('{0:s}'.format('%block ChemicalSpeciesLabel') + "\n")

    for i in range(0, nspecies):
        fout_siesta.write('{0:d} {1:d} {2:s}'.format(i+1,  atomic_num[i],  atms_name[i])+ "\n")
        ntype_order[i] = i+1

    fout_siesta.write('{0:s}'.format('%endblock ChemicalSpeciesLabel') + "\n")
    fout_siesta.write('{0:s}'.format('%block LatticeVectors') + "\n")

    for i in range(3):
        fout_siesta.write('{0:12.8f} {1:12.8f} {2:12.8f}'.format(latt_vec[i][0], latt_vec[i][1], latt_vec[i][2]) + "\n")
    
    fout_siesta.write('{0:s}'.format('%endblock LatticeVectors') + "\n")

    fout_siesta.write('{0:s}'.format('AtomicCoordinatesFormat                Fractional') + "\n")
    fout_siesta.write('{0:s}'.format('%block AtomicCoordinatesAndAtomicSpecies') + "\n")

    ncount = 0

    ntype = 0
    for i in range(0, natms[ntype]):
        ncount = ncount + 1

        m = int(bi_atms[i][0])
        
        x = bi_atms[i][1]
        y = bi_atms[i][2]
        z = bi_atms[i][3]
        fout_siesta.write('{:12.8f} {:12.8f} {:12.8f} {:d} {:d} '.format(x/a, y/b, z/c, ntype_order[ntype], atomic_num[ntype]) +   "\n")
       
    ntype = 1
    
    fout_siesta_spin.write('{0:s}'.format('%block DM.InitSpin')+"\n")

    for i in range(0, natms[ntype]):

        ncount = ncount + 1    
        m = int(fe_atms[i][0])
        
        x = fe_atms[i][1]
        y = fe_atms[i][2]
        z = fe_atms[i][3]
        if m ==1:
          fout_siesta.write('{:12.8f} {:12.8f} {:12.8f} {:d} {:d} '.format(x/a, y/b, z/c, ntype_order[ntype], atomic_num[ntype]) +   "\n")
          fout_siesta_spin.write('{0:d} {1:s}'.format(ncount,  '  +')+"\n")
        else:
          fout_siesta.write('{:12.8f} {:12.8f} {:12.8f} {:d} {:d} '.format(x/a, y/b, z/c, ntype_order[ntype], atomic_num[ntype]) +   "\n")
          fout_siesta_spin.write('{0:d} {1:s}'.format(ncount,  '  -')+"\n")

    fout_siesta_spin.write('{0:s}'.format('%endblock DM.InitSpin')+"\n")

    ntype = 2
  
    for i in range(0, natms[ntype]):
        ncount = ncount + 1

        m = int(o_atms[i][0])
        
        x = o_atms[i][1]
        y = o_atms[i][2]
        z = o_atms[i][3]
        fout_siesta.write('{:12.8f} {:12.8f} {:12.8f} {:d} {:d} '.format(x/a, y/b, z/c, ntype_order[ntype], atomic_num[ntype]) +   "\n")
                
    fout_siesta.write('{0:s}'.format('%endblock AtomicCoordinatesAndAtomicSpecies') + "\n")

#
def qe_output(natms, latt_vec, bi_atms,  fe_atms, o_atms):

    a = latt_vec[0][0] # [1-10]
    b = latt_vec[1][1] # [110]
    c = latt_vec[2][2] # [001]

    fout_qe = open('qe.in', 'w')
    fout_qe.write('CELL_PARAMETERS angstrom' + "\n")

    for i in range(3):
        fout_qe.write('{0:12.8f} {1:12.8f} {2:12.8f}'.format(latt_vec[i][0], latt_vec[i][1], latt_vec[i][2]) + "\n")
    fout_qe.write("\n\n" + 'ATOMIC_POSITIONS (crystal)' + "\n")

    ncount = 0

    ntype = 0
    for i in range(0, natms[ntype]):
        ncount = ncount + 1

        m = int(bi_atms[i][0])
        
        x = bi_atms[i][1]
        y = bi_atms[i][2]
        z = bi_atms[i][3]
        dispx = 1
        dispy = 1
        dispz = 1
        fout_qe.write('Bi  ' + '{:12.8f} {:12.8f} {:12.8f} {:d} {:d} {:d}'.format(x/a, y/b, z/c, dispx, dispy, dispz) +   "\n")
       
    ntype = 1
    for i in range(0, natms[ntype]):

        ncount = ncount + 1    
        m = int(fe_atms[i][0])
        
        x = fe_atms[i][1]
        y = fe_atms[i][2]
        z = fe_atms[i][3]
        dispx = 1
        dispy = 1
        dispz = 1

        if m ==1:
           fout_qe.write('Fe1  ' + '{:12.8f} {:12.8f} {:12.8f} {:d} {:d} {:d}'.format(x/a, y/b, z/c, dispx, dispy, dispz) +   "\n")
        else:
           fout_qe.write('Fe2  ' + '{:12.8f} {:12.8f} {:12.8f} {:d} {:d} {:d}'.format(x/a, y/b, z/c, dispx, dispy, dispz) +   "\n")	
 
   
    ntype = 2
    for i in range(0, natms[ntype]):
        ncount = ncount + 1

        m = int(o_atms[i][0])
        
        x = o_atms[i][1]
        y = o_atms[i][2]
        z = o_atms[i][3]
        dispx = 1
        dispy = 1
        dispz = 1
        fout_qe.write('O  ' + '{:12.8f} {:12.8f} {:12.8f} {:d} {:d} {:d}'.format(x/a, y/b, z/c, dispx, dispy, dispz) +   "\n")
      

#
def vasp_output(nspecies, natms, atms_name, latt_vec, bi_atms,  fe_atms, o_atms):
    fout_vasp = open("input_struc.vasp", 'w')
    # print(latt_vec)
    fout_vasp.write('Cs Pb Cl' + "\n")
    fout_vasp.write("1.0" + "\n")
    for i in range(3):
        # print(latt_vec[i][0])
        fout_vasp.write('{:12.8f} {:12.8f} {:12.8f}'.format(latt_vec[i][0], latt_vec[i][1], latt_vec[i][2]) + "\n")
    

    for i in range(0, nspecies):
        fout_vasp.write('{:s}'.format(atms_name[i])+ '  ')   
    fout_vasp.write("\n")
    for i in range(0, nspecies):
        fout_vasp.write('{:d}'.format(natms[i]) + '  ')   
    fout_vasp.write("\n")
    fout_vasp.write("Cartesian" + "\n")

    ntype = 0
    for i in range(0, natms[ntype]):

        #m = int(bi_atms[i][0])
        
        x = bi_atms[i][0]
        y = bi_atms[i][1]
        z = bi_atms[i][2]
        fout_vasp.write('{0:12.8f} {1:12.8f} {2:12.8f}'.format(x, y, z) + "\n")
       
    ntype = 1
   
    for i in range(0, natms[ntype]):
        #m = int(fe_atms[i][0])
        
        x = fe_atms[i][0]
        y = fe_atms[i][1]
        z = fe_atms[i][2]
#        if m ==1:
        fout_vasp.write('{0:12.8f} {1:12.8f} {2:12.8f}'.format(x, y, z) + "\n")
 #       else:
 #        fout_vasp.write('{0:12.8f} {1:12.8f} {2:12.8f}'.format(x, y, z) + "\n")

    ntype = 2
   
    for i in range(0, natms[ntype]):

        #m = int(o_atms[i][0])
        
        x = o_atms[i][0]
        y = o_atms[i][1]
        z = o_atms[i][2]
        fout_vasp.write('{0:12.8f} {1:12.8f} {2:12.8f}'.format(x, y, z) + "\n")
                
