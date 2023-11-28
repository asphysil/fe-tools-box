import numpy as np 

from Ctype_AFM_order import *
from nearest_neighbour_GtypeAFM import *

from nearest_neighbour_2nd import *

# getting required information for in-phase o6 rotation 
def inphase_info(axis, nx, ny, nz, nubi, nufe, nuo, latt_unit, bi_atmu, fe_atmu, o_atmu):
    # C-type AFM order 
    if axis==1: # x-axis 
        fe_ctypeAFM=CtypeAFM_order_xaxis(nx, ny, nz, nufe, latt_unit, fe_atmu) #along x-axis and yz-plane AFM
    elif axis==2: # y-axis
        fe_ctypeAFM=CtypeAFM_order_yaxis(nx, ny, nz, nufe, latt_unit, fe_atmu) # along y-axis and xz-plane AFM
    else: #z-axis
        fe_ctypeAFM=CtypeAFM_order_zaxis(nx, ny, nz, nufe, latt_unit, fe_atmu)  # along z-axis and xy-plane AFM
    #print("-----")
    #print("-----")
    # ordering and arranging atoms 
    # 4 columns data 
    fe_ctypeAFM[:, [4, 1]] = fe_ctypeAFM[:, [1, 4]] # Column swapped 
    fe_ctypeAFM=fe_ctypeAFM[fe_ctypeAFM[:,4].argsort()] 
    fe_ctypeAFM[:, [4, 1]] = fe_ctypeAFM[:, [1, 4]] # Column swapped 
    fe_ctypeAFM=fe_ctypeAFM[:, 1:5] # Ctype-AFM
    #print(fe_ctypeAFM)
    
    nn = 6  # Each B atoms has 6 Oxygen atoms
    dbo = 4.0/2.0 + 1.5
    print("********Warning******")
    print("**Please change if B-O bond length is not equal to ", dbo, '  for you system **') 
    #atms_bo= BO_1stdist_pbc(nn, nufe, nuo, latt_unit, dbo,  fe_atms_order, o_atmu)
    
    atms_bo= BO_1stdist_nonpbc(nn, nufe, nuo, dbo, fe_ctypeAFM, o_atmu) 
    #print('------', atms_bo)
    #exit()
    nn_ab = 8
    dab = 4.0*(np.sqrt(3)/2.0) + 1.0 
    print("********Warning******")
    print("**Please change if A-B bond length is not equal to ", dab, '  for you system **') 
    atms_ab, coord_ab = AB_1stdist_pbc(nn_ab, nubi, nufe, dab, latt_unit, bi_atmu, fe_ctypeAFM)
    #print('-----------------------------------')
    #print(atms_ab)
    #print(coord_ab)
    local_axis = cal_local_axis(nn_ab, nufe, coord_ab)
    
    nn_bo=6
    loc_noatm, loc_axis_ordered = local_axis_arrange(nn_bo, nufe, local_axis, atms_bo)
    return (loc_axis_ordered, loc_noatm, atms_bo)


def antiphase_info(nx, ny, nz, nubi, nufe, nuo, latt_unit, bi_atmu, fe_atmu, o_atmu):

    nn = 6 # Each Fe atoms has 6 Fe atoms
    fe_GtypeAFM=GtypeAFE_order_v2(nn, nufe, latt_unit, fe_atmu)
    ###############################################
    #print("----")
    # ordering and arranging atoms 
    # 3  columns data 
    fe_GtypeAFM[:, [3, 0]] = fe_GtypeAFM[:, [0, 3]] # Column swapped 0 and z
    fe_GtypeAFM=fe_GtypeAFM[fe_GtypeAFM[:,3].argsort()] 
    fe_GtypeAFM[:, [3, 0]] = fe_GtypeAFM[:, [0, 3]] # Column swapped 0 and z
    #print(fe_atms_order_g)
    
    # G-type AFM Order
    nn = 6
    dbo = 4.0/2.0 + 1.5 
    print("********Warning******")
    print("**Please change if B-O bond length is not equal to ", dbo, '  for you system **') 
    atms_bo= BO_1stdist_nonpbc(nn, nufe, nuo, dbo, fe_GtypeAFM, o_atmu) 
    #print('------', atms_bo)
    #exit()
    nn_ab = 8
    dab = 4.0*(np.sqrt(3)/2.0) + 1.0 
    print("********Warning******")
    print("**Please change if A-B bond length is not equal to ", dab, '  for you system **') 
    atms_ab, coord_ab= AB_1stdist_pbc(nn_ab, nubi, nufe, dab, latt_unit, bi_atmu, fe_GtypeAFM)
   #print('-----------------------------------')
    #print(atms_ab)
    #print(coord_ab)
    local_axis= cal_local_axis(nn_ab, nufe, coord_ab)
    
    nn_bo=6
    loc_noatm, loc_axis_ordered= local_axis_arrange(nn_bo, nufe, local_axis, atms_bo)
    return  (loc_axis_ordered, loc_noatm, atms_bo)