import numpy as np

from O6_rotation import *
from SubFunction import * 

from lm_inputs import *
from lm_outputs import * 

from nearest_neighbour_2nd import *


from single_domain import *
from lammps_output import * 

#from Ctype_AFM_order import *

# reading input
file_name ='POSCAR'

latt_unit, natms, bi_atmu, fe_atmu, o_atmu = read_input_vasp_ABO3(file_name)

a = latt_unit[0][0] # [1-10]
b = latt_unit[1][1] # [110]
c = latt_unit[2][2] # [001]

nubi = natms[0]
nufe = natms[1]
nuo = natms[2]

# Supper cell dimension of a cubic unit cell of ABO3
print("Please check supper cell dimension")
nx=2
ny=2
nz=2
#
if (nx*ny*nz !=nufe):
    print("Supercell dimension is wrong, please change it")
    exit()

#
#fe=CtypeAFM_order_xaxis(nx, ny, nz, nufe, latt_unit, fe_atmu)
#fe=CtypeAFM_order_yaxis(nx, ny, nz, nufe, latt_unit, fe_atmu)
print("----")

#exit()

nn=6
dbo = 4.0/2.0 + 1.5

phase=[-1, 1,-1] # -1=antiphase and 1=inphase
rot_angle=[5,5,5] # Rotation Angle

o_atms_rot=np.copy(o_atmu)  # copy the intial values

for i in range(1,4): # Rotation axis 
    axis=i
    if phase[i-1]==1:
        loc_axis_ordered, loc_noatm, atms_bo=inphase_info(axis, nx, ny, nz, nubi, nufe, nuo, latt_unit, bi_atmu, fe_atmu, o_atmu)
        alpha=rot_angle[i]
        o_atms_rot=O6_rotation_ABO3(nuo, nufe, axis, alpha, loc_axis_ordered, loc_noatm, atms_bo, o_atms_rot)
    else:
        loc_axis_ordered, loc_noatm, atms_bo=antiphase_info(nx, ny, nz, nubi, nufe, nuo, latt_unit, bi_atmu, fe_atmu, o_atmu)
        alpha=rot_angle[i-1]
        o_atms_rot=O6_rotation_ABO3(nuo, nufe, axis, alpha, loc_axis_ordered, loc_noatm, atms_bo, o_atms_rot)

    o_atms_rot=np.copy(o_atms_rot) # Update the calculated value



print("-----org---")
print(o_atmu)
print("---------")
print(o_atms_rot)


print("-------- No O atoms rotation")
center_of_mass_O6(latt_unit, nn, nufe, nuo, dbo, fe_atmu, o_atmu)	
print("--------")
print("--------With O atoms rotation")
center_of_mass_O6(latt_unit, nn, nufe, nuo, dbo, fe_atmu, o_atms_rot)	
print("--------")


fe_GtypeAFM=GtypeAFE_order_v2(nn, nufe, latt_unit, fe_atmu)
fe_GtypeAFM[:, [3, 0]] = fe_GtypeAFM[:, [0, 3]] # Column swapped 0 and z
fe_GtypeAFM=fe_GtypeAFM[fe_GtypeAFM[:,3].argsort()] 
fe_GtypeAFM[:, [3, 0]] = fe_GtypeAFM[:, [0, 3]] # Column swapped 0 and z

#exit()

#
#-----initial domain---------------
#disp = 0.22
#fe_atms_sd = sd_rombohedral(nufe,  disp, fe_atms_order)
#fe_atms_order = np.copy(fe_atms_sd)

disp = 0.52
#print('dx=', disp/np.sqrt(3))
bi_atms_sd = sd_rombohedral(nubi, disp, latt_unit, bi_atmu)



# #
# disp = -0.2
# fe_atms_sd = sd_tetragonal(nufe,  disp, fe_atms_order)
# fe_atms_order = np.copy(fe_atms_sd)

# disp = -0.51
# bi_atms_sd = sd_tetragonal(nubi, disp,  bi_atms_order)
# bi_atms_order = np.copy(bi_atms_sd)
#-----------------------------------

print("--------With O atoms rotation and Bi atoms displace")
center_of_mass_O12(latt_unit, nn, nubi, nuo, dbo, bi_atmu, o_atms_rot)

nspecies =3 
atms_name = ['Cs', 'Pb', 'Cl']

#print(fe_atmu)
vasp_output(nspecies, natms, atms_name, latt_unit, bi_atms_sd,  fe_atmu, o_atms_rot)
lammps_bv(natms, latt_unit, bi_atmu, fe_GtypeAFM, o_atms_rot)
