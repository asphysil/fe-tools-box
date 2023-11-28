import numpy as np



def read_input_vasp_ABO3(file_name):

    fread=open(file_name, 'r')
    #flines = fread.readlines()
    l=0
    i=0
    latt_unit = np.zeros((3,3))
    
    for lines in fread.readlines():
        l = l+1
	#print(lines)
        if l>2 and l<6:
           latt_unit[i] = [float(x) for x in lines.split()]
           i = i +1
        if l==6:
           print('atom order in .vasp file', lines)
        if l==7:
           ntype = [int(x) for x in lines.split()]
        if l==8:
           csys =lines[0].split()

    n = len(ntype)
    natms = np.zeros(n, dtype=int)
    for i in range(0, n):
        natms[i] = ntype[i]

    if csys[0]=='C' or csys=='c':
       print('Coorect coordinate system "Cartesian"')
       all_atms = np.genfromtxt(file_name, skip_header=8)

       n1 = int(ntype[0])
       n2 = int(ntype[1])
       n3= int(ntype[2])
  
       print('n1=', n1, 'n2=', n2, 'n3=', n3)
       bi_atmu   =all_atms[0:n1]
       m1 = n1+n2
       fe_atmu  =all_atms[n1:m1]
       m2 = m1+n3
       o_atmu =all_atms[m1:m2]
       
    else:
       print('Coordinate system is Direct or crystal')
       all_atms = np.genfromtxt(file_name, skip_header=8)

       n1 = int(ntype[0])
       n2 = int(ntype[1])
       n3= int(ntype[2])

       print('n1=', n1, 'n2=', n2, 'n3=', n3 )
       bi_atmu   =all_atms[0:n1]
       m1 = n1+n2
       fe_atmu  =all_atms[n1:m1]
       m2 = m1+n3
       o_atmu =all_atms[m1:m2]
      

       
       o_atmu  = np.matmul(np.copy(o_atmu), latt_unit)  # crystal to cartesian
       bi_atmu  = np.matmul(np.copy(bi_atmu), latt_unit)  # crystal to cartesian
       fe_atmu  = np.matmul(np.copy(fe_atmu), latt_unit)  # crystal to cartesian

    natms[0]=n1# Bi
    natms[1]=n2 #fe
    natms[2]=n3 # O
    return (latt_unit, natms,bi_atmu, fe_atmu, o_atmu)