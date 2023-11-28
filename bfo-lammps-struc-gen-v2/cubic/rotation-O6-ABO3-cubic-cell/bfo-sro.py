import numpy as np
from ase.io import read,write


from ase.build import cut,stack

import spglib as spg


    
instr1 = read("bfo.vasp", format="vasp")




#substrate_str1 =  read("sro.vasp", format="vasp")

#cell_sto = substrate_str1.get_cell()

# instr1.set_cell(cell_sto)
nlbfo=4
# nl=6
outstr1 = cut(instr1, a=(1,-1,0), b=(1,1,0), nlayers=nlbfo,)

# outstr2 = cut(substrate_str1, a=(1,-1,0), b=(1,1,0), nlayers=nl)

# interface = stack(outstr2, outstr1, maxstrain=1, distance=2.5)
write("str_out.cif", outstr1)
#fout='str-out-bfo-sro-' + str(nlbfo) + '-' + str(nl) + 'L.cif'
#write(fout, interface)

print(spg.get_version())
print(spg.get_spacegroup(outstr1))

