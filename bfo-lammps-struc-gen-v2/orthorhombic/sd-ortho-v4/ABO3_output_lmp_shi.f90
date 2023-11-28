
SUBROUTINE lammps_struc_inp_liu(nsba, nsti, nso1, nso2, nso3)


  USE my_constants 
  USE data_structure, only: latt_sup, a_atms, b_atms, o1_atms, o2_atms, o3_atms
                      
 IMPLICIT NONE 

INTEGER(dp_int), INTENT(IN) :: nsba, nsti, nso1, nso2, nso3

REAL(dp_real) :: sqr2
REAL(dp_real), dimension(4, nsti):: a1_atms, b1_atms, o11_atms, o22_atms, o33_atms
REAL (dp_real) :: latt_strain(3,3)

INTEGER(dp_int) :: ntot_atms,  ncout
INTEGER(dp_int) :: i, j


ntot_atms = nsba + nsti + nso1 + nso2 + nso3



DO i =1, 3
DO j =1, 3
latt_strain(i,j) =0.0
ENDDO
ENDDO

sqr2 = SQRT(2.0)

latt_strain(1,1) = (a11 + (strainx*a11*0.01))*DBLE(nx)
latt_strain(2,2) = (a22 + (strainy*a22*0.01))*DBLE(ny)
latt_strain(3,3) = (a33 + (strainz*a33*0.01))*DBLE(nz)

!DO i =1, 3
!    PRINT*, "Latt s=0", (latt_sup(i,j), j=1,3)
!ENDDO

!DO i =1, 3
!    PRINT*, "Latt s/=0", (latt_strain(i,j), j=1,3)
!ENDDO

DO i =1, nsti
   a1_atms(1,i) = a_atms(1,i)
   b1_atms(1,i) = b_atms(1,i)
  o11_atms(1,i) =o1_atms(1,i)
  o22_atms(1,i) =o2_atms(1,i)
  o33_atms(1,i) =o3_atms(1,i)
  
     DO j = 1,3
      a1_atms(j+1,i) =  (a_atms(j+1,i)/latt_sup(j,j))*latt_strain(j,j)
      b1_atms(j+1,i) =  (b_atms(j+1,i)/latt_sup(j,j))*latt_strain(j,j)
     o11_atms(j+1,i) = (o1_atms(j+1,i)/latt_sup(j,j))*latt_strain(j,j)
     o22_atms(j+1,i) = (o2_atms(j+1,i)/latt_sup(j,j))*latt_strain(j,j)
     o33_atms(j+1,i) = (o3_atms(j+1,i)/latt_sup(j,j))*latt_strain(j,j)
     ENDDO

ENDDO

!!!!!!!!!!
! domain
!SD
CALL domain_struc_A(nsba, a1_atms)
CALL domain_struc_B(nsti, b1_atms) 
!PRINT*, 'ok'

OPEN(UNIT=15, FILE='data-shi.BTO', ACTION='WRITE')    
    !!!!!!!!!!!!!!!!!!!! Writing lammps dada file format!!!!!!!!!
WRITE(15,*)'#LAMMPS ABO3 domain'
WRITE(15,*)"  "
WRITE(15,190) ntot_atms, "atoms"
!WRITE(15,191) ntot_angs, "angles"

!PRINT*, " Supercell(nx*ny*nz)=", nx, ny, nz
!PRINT*, " Total number of atoms=", ntot_atms
!PRINT*,"**************Warning:********************"
!PRINT*, "If total number of atoms greater than ten million atoms, Please modify FORMAT and INTEGER options" 
!PRINT*, "                               "

190 FORMAT(I9, A6)
!191 FORMAT(I9, 1X, A6)

WRITE(15,192) 3, "atom", "types"
!WRITE(15,192) 1, "angle", "types"
WRITE(15,*)"  "

192 FORMAT(I2, A6, A6)

WRITE(15,200)0.0, latt_strain(1,1),  "xlo xhi"
WRITE(15,200)0.0, latt_strain(2,2),  "ylo yhi"
WRITE(15,200)0.0, latt_strain(3,3),  "zlo zhi"
200   FORMAT(1X, F3.1, F12.5, 2X, A7)

WRITE(15,*)"  "
WRITE(15,*)"Masses"
WRITE(15,*)"  "

!PRINT*,"**************Warning:********************"
!PRINT*, "if the structure is not BaTiO3 Please modify atomic masses" 
!PRINT*, "                                             "

WRITE(15,202)"1", massA ! " # Ba" !Modify this for other ABO3 compounds
WRITE(15,202)"2", massB !" # Ti"   ! Modify this for other ABO3 compounds
WRITE(15,202)"3", massO !" # O" ! Modify this for other ABO3 compounds
WRITE(15,202)"   "
!201 FORMAT(I2, F9.4)
202 FORMAT(A2, F12.5)

!WRITE(15,*)"Angle Coeffs"
!WRITE(15,*)"  "

!PRINT*,"**************Warning:********************"
!PRINT*, "if the structure is not BaTiO3 Please modify Angle Coeffs "
!PRINT*, "                                                        "

!WRITE(15,205)1, angle_coeffs,  180.00 ! Modify this for other ABO3 compounds
!205 FORMAT(1X,I1, 1X, I3, 1X, F7.2)

WRITE(15,*)"  "
WRITE(15,*)"Atoms"
WRITE(15,*)"  "
!Modify this for other ABO3 compounds

!PRINT*,"**************Warning:********************"
!PRINT*, "if the structure is not BaTiO3 Please modify atomic charges"
!PRINT*, "                                                        "

!charge_A=1.34730 
!charge_B=1.28905
!charge_O=-0.87878
!

ncout = 0
  DO i=1, nsba 
   ncout = ncout + 1
   WRITE(15,100)ncout, 1, 1, charge_A,(a1_atms(j, i),j=2,4)
   !PRINT*, "ok"
    ENDDO


DO i=1, nsti
!  PRINT*, "ok"
  ncout = ncout + 1
   WRITE(15,100)ncout, 1, 2, charge_B, (b1_atms(j,i),j=2,4)
ENDDO


DO i=1, nso1 
  ncout = ncout + 1
  WRITE(15,100)ncout, 1, 3, charge_O, (o11_atms(j,i),j=2,4)
ENDDO

DO i=1, nso2
  ncout = ncout + 1
  WRITE(15,100)ncout, 1, 3, charge_O,  (o22_atms(j,i),j=2,4)
  ENDDO

DO i=1, nso3
  ncout = ncout + 1
  WRITE(15,100)ncout, 1, 3, charge_O, (o33_atms(j,i),j=2,4)
  ENDDO

100 FORMAT(I10, 1X, I3, 1X, I3, 1X, F10.5, 1X, 3F12.7)
END SUBROUTINE lammps_struc_inp_liu


SUBROUTINE lammps_struc_inp_liu_2nd(nsba, nsti, nso1, nso2, nso3)


  USE my_constants 
  USE data_structure, only: latt_sup, a_atms, b_atms, o_atms
                      
 IMPLICIT NONE 

INTEGER(dp_int), INTENT(IN) :: nsba, nsti, nso1, nso2, nso3

REAL(dp_real) :: sqr2
REAL(dp_real), dimension(4, nsti):: a1_atms, b1_atms

REAL(dp_real), dimension(4, 3*nsti):: o11_atms

REAL (dp_real) :: latt_strain(3,3)

INTEGER(dp_int) :: ntot_atms,  ncout, nso 
INTEGER(dp_int) :: i, j


ntot_atms = nsba + nsti + nso1 + nso2 + nso3

nso = nso1 + nso2 + nso3

DO i =1, 3
DO j =1, 3
latt_strain(i,j) =0.0
ENDDO
ENDDO

sqr2 = SQRT(2.0)

latt_strain(1,1) = (a11 + (strainx*a11*0.01))*DBLE(nx)
latt_strain(2,2) = (a22 + (strainy*a22*0.01))*DBLE(ny)
latt_strain(3,3) = (a33 + (strainz*a33*0.01))*DBLE(nz)

!DO i =1, 3
!    PRINT*, "Latt s=0", (latt_sup(i,j), j=1,3)
!ENDDO

!DO i =1, 3
!    PRINT*, "Latt s/=0", (latt_strain(i,j), j=1,3)
!ENDDO

DO i =1, nsti
   a1_atms(1,i) = a_atms(1,i)
   b1_atms(1,i) = b_atms(1,i)
     DO j = 1,3
      a1_atms(j+1,i) =  (a_atms(j+1,i)/latt_sup(j,j))*latt_strain(j,j)
      b1_atms(j+1,i) =  (b_atms(j+1,i)/latt_sup(j,j))*latt_strain(j,j)
     ENDDO

ENDDO

DO i =1 , nso 
  o11_atms(1,i) = o_atms(1,i)
  DO j =1 , 3
    o11_atms(j+1,i) = (o_atms(j+1,i)/latt_sup(j,j))*latt_strain(j,j)
  ENDDO 
ENDDO 
!!!!!!!!!!
! domain
!SD
CALL domain_struc_A(nsba, a1_atms)
CALL domain_struc_B(nsti, b1_atms) 
!PRINT*, 'ok'

OPEN(UNIT=15, FILE='data-shi.BTO', ACTION='WRITE')    
    !!!!!!!!!!!!!!!!!!!! Writing lammps dada file format!!!!!!!!!
WRITE(15,*)'#LAMMPS ABO3 domain'
WRITE(15,*)"  "
WRITE(15,190) ntot_atms, "atoms"
!WRITE(15,191) ntot_angs, "angles"

!PRINT*, " Supercell(nx*ny*nz)=", nx, ny, nz
!PRINT*, " Total number of atoms=", ntot_atms
!PRINT*,"**************Warning:********************"
!PRINT*, "If total number of atoms greater than ten million atoms, Please modify FORMAT and INTEGER options" 
!PRINT*, "                               "

190 FORMAT(I9, A6)
!191 FORMAT(I9, 1X, A6)

WRITE(15,192) 3, "atom", "types"
!WRITE(15,192) 1, "angle", "types"
WRITE(15,*)"  "

192 FORMAT(I2, A6, A6)

WRITE(15,200)0.0, latt_strain(1,1),  "xlo xhi"
WRITE(15,200)0.0, latt_strain(2,2),  "ylo yhi"
WRITE(15,200)0.0, latt_strain(3,3),  "zlo zhi"
200   FORMAT(1X, F3.1, F12.5, 2X, A7)

WRITE(15,*)"  "
WRITE(15,*)"Masses"
WRITE(15,*)"  "

!PRINT*,"**************Warning:********************"
!PRINT*, "if the structure is not BaTiO3 Please modify atomic masses" 
!PRINT*, "                                             "

WRITE(15,202)"1", massA ! " # Ba" !Modify this for other ABO3 compounds
WRITE(15,202)"2", massB !" # Ti"   ! Modify this for other ABO3 compounds
WRITE(15,202)"3", massO !" # O" ! Modify this for other ABO3 compounds
WRITE(15,202)"   "
!201 FORMAT(I2, F9.4)
202 FORMAT(A2, F12.5)

!WRITE(15,*)"Angle Coeffs"
!WRITE(15,*)"  "

!PRINT*,"**************Warning:********************"
!PRINT*, "if the structure is not BaTiO3 Please modify Angle Coeffs "
!PRINT*, "                                                        "

!WRITE(15,205)1, angle_coeffs,  180.00 ! Modify this for other ABO3 compounds
!205 FORMAT(1X,I1, 1X, I3, 1X, F7.2)

WRITE(15,*)"  "
WRITE(15,*)"Atoms"
WRITE(15,*)"  "
!Modify this for other ABO3 compounds

!PRINT*,"**************Warning:********************"
!PRINT*, "if the structure is not BaTiO3 Please modify atomic charges"
!PRINT*, "                                                        "

!charge_A=1.34730 
!charge_B=1.28905
!charge_O=-0.87878
!

ncout = 0
  DO i=1, nsba 
   ncout = ncout + 1
   WRITE(15,100)ncout, 1, 1, charge_A,(a1_atms(j, i),j=2,4)
   !PRINT*, "ok"
    ENDDO


DO i=1, nsti
!  PRINT*, "ok"
  ncout = ncout + 1
   WRITE(15,100)ncout, 1, 2, charge_B, (b1_atms(j,i),j=2,4)
ENDDO


DO i=1, nso
  ncout = ncout + 1
  WRITE(15,100)ncout, 1, 3, charge_O, (o11_atms(j,i),j=2,4)
ENDDO


100 FORMAT(I10, 1X, I3, 1X, I3, 1X, F10.5, 1X, 3F12.7)
END SUBROUTINE lammps_struc_inp_liu_2nd