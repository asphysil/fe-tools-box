
SUBROUTINE struc_vasp_format(nsba, nsti, nso1, nso2, nso3)

  USE my_constants 
  USE data_structure, only: latt_sup, a_atms, b_atms, o1_atms, o2_atms, o3_atms
 IMPLICIT NONE 

INTEGER(dp_int), INTENT(IN) :: nsba, nsti, nso1, nso2, nso3


REAL(dp_real) :: sqr2
REAL(dp_real), dimension(4, nsti):: a1_atms, b1_atms, o11_atms, o22_atms, o33_atms
REAL (dp_real) :: latt_strain(3,3)
    


INTEGER(dp_int) :: i, j
    

sqr2= SQRT(2.0)

!!!!!!!!!
DO i =1, 3
DO j =1, 3
latt_strain(i,j) =0.0
ENDDO
ENDDO


latt_strain(1,1) = (a11 + (strainx*a11*0.01))*DBLE(nx)
latt_strain(2,2) = (a22 + (strainy*a22*0.01))*DBLE(ny)
latt_strain(3,3) = (a33 + (strainz*a33*0.01))*DBLE(nz)

!PRINT*, latt_strain(1,1), latt_strain(2,2), latt_strain(3,3), (latt_sup(j,j), j=1,3)

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

!

OPEN(UNIT=15, FILE='struc_out.vasp', ACTION='WRITE')    
    !!!!!!!!!!!!!!!!!!!! Writing lammps dada file format!!!!!!!!!
WRITE(15,*)'# ABO3 domain'
WRITE(15,*)" 1.0"

!PRINT*, " Supercell(nx*ny*nz)=", nx, ny, nz
!PRINT*, " Total number of atoms=", (nsba + nsti + nso1 + nso2 + nso3)
!PRINT*,"**************Warning:********************"
!PRINT*, "If total number of atoms greater than ten million atoms, Please modify FORMAT and INTEGER options" 
!PRINT*, "                               "

DO i = 1, 3
WRITE(15,200)latt_strain(i,1), latt_strain(i,2), latt_strain(i,3)
ENDDO
200 FORMAT(1X, 3F12.8)
!200   FORMAT(1X, F3.1, F12.5, 2X, A7)

WRITE(15,*)" Bi   Fe   O "

WRITE(15,201)nsba, nsti, nso1+nso2+nso3
201 FORMAT(1X, 3I5) 

!PRINT*,"**************Warning:********************"
!PRINT*, "if the structure is not BaTiO3 Please modify atomic masses" 
!PRINT*, "                                             "

!PRINT*,"**************Warning:********************"
!PRINT*, "if the structure is not BaTiO3 Please modify Angle Coeffs "
!PRINT*, "                                                        "


WRITE(15,*)"Cartesian"
!Modify this for other ABO3 compounds
!PRINT*,"**************Warning:********************"
!PRINT*, "if the structure is not BaTiO3 Please modify atomic charges"
!PRINT*, "                                                        "

  DO i=1, nsba 
   WRITE(15,100)(a1_atms(j,i),j=2,4)
   !PRINT*, i, (a1_atms(i, j),j=2,4)
    ENDDO


DO i=1, nsti
!  PRINT*, "ok"
   WRITE(15,100)(b1_atms(j,i),j=2,4)
ENDDO

DO i=1, nso1 
  WRITE(15,100) (o11_atms(j,i),j=2,4)
ENDDO

DO i=1, nso2
  WRITE(15,100)(o22_atms(j,i),j=2,4)
  ENDDO

DO i=1, nso3
  WRITE(15,100)(o33_atms(j,i),j=2,4)
  ENDDO

100 FORMAT(1X, 3F12.7)
END SUBROUTINE struc_vasp_format


SUBROUTINE struc_vasp_format_2nd(nsba, nsti, nso1, nso2, nso3)

  USE my_constants 
  USE data_structure, only: latt_sup, a_atms, b_atms, o_atms
 IMPLICIT NONE 

INTEGER(dp_int), INTENT(IN) :: nsba, nsti, nso1, nso2, nso3


REAL(dp_real) :: sqr2
REAL(dp_real), dimension(4, nsti):: a1_atms, b1_atms
REAL(dp_real), dimension(4, 3*nsti):: o11_atms
REAL (dp_real) :: latt_strain(3,3)
    


INTEGER(dp_int) :: i, j, nso 
    

nso = nso1 + nso2 + nso3 

sqr2= SQRT(2.0)

!!!!!!!!!
DO i =1, 3
DO j =1, 3
latt_strain(i,j) =0.0
ENDDO
ENDDO


latt_strain(1,1) = (a11 + (strainx*a11*0.01))*DBLE(nx)
latt_strain(2,2) = (a22 + (strainy*a22*0.01))*DBLE(ny)
latt_strain(3,3) = (a33 + (strainz*a33*0.01))*DBLE(nz)

!PRINT*, latt_strain(1,1), latt_strain(2,2), latt_strain(3,3), (latt_sup(j,j), j=1,3)

DO i =1, nsti
  a1_atms(1,i) = a_atms(1,i)
  b1_atms(1,i) = b_atms(1,i)
  
     DO j = 1,3
      a1_atms(j+1,i) =  (a_atms(j+1,i)/latt_sup(j,j))*latt_strain(j,j)
      b1_atms(j+1,i) =  (b_atms(j+1,i)/latt_sup(j,j))*latt_strain(j,j)
ENDDO

ENDDO

DO i=1, nso 
  o11_atms(1,i) =o_atms(1,i)
  DO j =1, 3
    o11_atms(j+1,i) = (o_atms(j+1,i)/latt_sup(j,j))*latt_strain(j,j)
  ENDDO 
ENDDO

!!!!!!!!!!
! domain
!SD
CALL domain_struc_A(nsba, a1_atms)
CALL domain_struc_B(nsti, b1_atms) 

!

OPEN(UNIT=15, FILE='struc_out.vasp', ACTION='WRITE')    
    !!!!!!!!!!!!!!!!!!!! Writing lammps dada file format!!!!!!!!!
WRITE(15,*)'# ABO3 domain'
WRITE(15,*)" 1.0"

!PRINT*, " Supercell(nx*ny*nz)=", nx, ny, nz
!PRINT*, " Total number of atoms=", (nsba + nsti + nso1 + nso2 + nso3)
!PRINT*,"**************Warning:********************"
!PRINT*, "If total number of atoms greater than ten million atoms, Please modify FORMAT and INTEGER options" 
!PRINT*, "                               "

DO i = 1, 3
WRITE(15,200)latt_strain(i,1), latt_strain(i,2), latt_strain(i,3)
ENDDO
200 FORMAT(1X, 3F12.8)
!200   FORMAT(1X, F3.1, F12.5, 2X, A7)

WRITE(15,*)" Bi   Fe   O "

WRITE(15,201)nsba, nsti, nso1+nso2+nso3
201 FORMAT(1X, 3I5) 

!PRINT*,"**************Warning:********************"
!PRINT*, "if the structure is not BaTiO3 Please modify atomic masses" 
!PRINT*, "                                             "

!PRINT*,"**************Warning:********************"
!PRINT*, "if the structure is not BaTiO3 Please modify Angle Coeffs "
!PRINT*, "                                                        "


WRITE(15,*)"Cartesian"
!Modify this for other ABO3 compounds
!PRINT*,"**************Warning:********************"
!PRINT*, "if the structure is not BaTiO3 Please modify atomic charges"
!PRINT*, "                                                        "

  DO i=1, nsba 
   WRITE(15,100)(a1_atms(j,i),j=2,4)
   !PRINT*, i, (a1_atms(i, j),j=2,4)
    ENDDO


DO i=1, nsti
!  PRINT*, "ok"
   WRITE(15,100)(b1_atms(j,i),j=2,4)
ENDDO

DO i=1, nso 
  WRITE(15,100) (o11_atms(j,i),j=2,4)
ENDDO


100 FORMAT(1X, 3F12.7)
END SUBROUTINE struc_vasp_format_2nd