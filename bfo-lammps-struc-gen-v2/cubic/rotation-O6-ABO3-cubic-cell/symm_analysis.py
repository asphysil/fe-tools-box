import numpy as np
from ase.io import read,write


from ase.build import cut,stack

import spglib as spg


#file_in ='pm-3m.vasp'   
file_in = 'input_struc.vasp' 
instr1 = read(file_in, format="vasp")

print(spg.get_version())
print(spg.get_spacegroup(instr1, symprec=0.4))

