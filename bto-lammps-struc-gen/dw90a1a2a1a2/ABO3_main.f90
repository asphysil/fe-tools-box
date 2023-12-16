!!!!!!!!!!!!!!!!!Main Program!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PROGRAM ABO3_domain_struc
      USE dtype
      IMPLICIT NONE
      !INTEGER, PARAMETER :: dp_real=SELECTED_REAL_KIND(16)
      REAL(dp_real), ALLOCATABLE :: A_atoms(:,:), B_atoms(:,:),  &
                           O1_atoms(:,:), O2_atoms(:,:), O3_atoms(:,:), &
                           A_s(:,:), B_s(:,:), &
                           O1_s(:,:), O2_s(:,:), O3_s(:,:)
     


      REAL(dp_real) :: Ba_atoms(1,3), Ti_atoms(1,3), O_atoms(3,3), &
                       a_unit(3,3), a_strain(3,3),  b(3), bf(3), T(3), angle, &
                       Ba_o(1,3), Ti_o(1,3), O_o(3,3), Ba_s(1,3), Ti_s(1,3), O_s(3,3), T_s(3)  !,&
!                       charge_A, charge_B, charge_O
 
      INTEGER(dp_int) :: nxy, nb,  ntot, no1, no2, no3
      INTEGER(dp_int) :: n1o, n2o, n3o
      INTEGER(dp_int) :: i, j, k, n, m1  ! integer for loop 
     
PRINT*, "               "
!PRINT*, " Warning:-        "
!
!PRINT*,"\o/"
!PRINT*," |"
!PRINT*,"/ \"
!PRINT*, "  "
!
!PRINT*, " Upto ten million atoms"
PRINT*, "                           "

!      PRINT*, " Enter dimension of the supercell nx*ny*nz"
!      READ*, nx, ny, nz
!nx=70;ny=70;nz=70
ntot=nx*ny*nz
nxy=nx*ny

nb=ntot+1
no1=2*ntot+1
no2=3*ntot+1
no3=4*ntot+1

ALLOCATE(A_atoms(ntot,3), B_atoms(nb:2*ntot,3), &
         O1_atoms(no1:3*ntot,3), O2_atoms(no2:4*ntot,3), O3_atoms(no3:5*ntot,3))

ALLOCATE(A_s(ntot,3),B_s(nb:2*ntot,3), &
        O1_s(no1:3*ntot,3), O2_s(no2:4*ntot,3), O3_s(no3:5*ntot,3))

DO i=1,3
  Do j=1,3
!  a_sup(i,j)=0.0
!  a_sup_T(i,j)=0.0
  a_unit(i,j)=0.0
  a_strain(i,j)=0.0
  ENDDO
    ENDDO

!      OPEN(UNIT=12, FILE='POSCAR', ACTION='READ')
      OPEN(UNIT=13, FILE='SUPERCELL.vasp', ACTION='WRITE')
      OPEN(UNIT=14, FILE='SUPERCELL_angles.dat', ACTION='WRITE')
      OPEN(UNIT=15, FILE='data.BTO', ACTION='WRITE')
PRINT*,       
!PRINT*, " Reading 5 atoms unit cell "
     ! read junck line!
!      READ(12,*)
!      READ(12,*) 
!      DO i =1, 3
!        READ(12,*) a_unit(i,1), a_unit(i,2),a_unit(i,3)
!        ENDDO
!      READ(12,*)
!      READ(12,*)
!      READ(12,*)
!      
!     READ(12,*) Ba_atoms(1,1), Ba_atoms(1,2), Ba_atoms(1,3)
!     READ(12,*) Ti_atoms(1,1),Ti_atoms(1,2),Ti_atoms(1,3)
!        DO i=1,3
!           READ(12,*) O_atoms(i,1),O_atoms(i,2),O_atoms(i,3)
!           ENDDO

a_unit(1,1)=a11
a_unit(2,2)=a22
a_unit(3,3)=a33

a_strain(1,1)=a_unit(1,1)+ (strainx*a_unit(1,1)*0.01)
a_strain(2,2)=a_unit(2,2)+ (strainy*a_unit(2,2)*0.01)
a_strain(3,3)=a_unit(3,3)+ (strainz*a_unit(3,3)*0.01)

Ba_atoms(1,1)=0.0; Ba_atoms(1,2)=0.0; Ba_atoms(1,3)=0.0  ! Ba atoms
Ti_atoms(1,1)=0.5; Ti_atoms(1,2)=0.5; Ti_atoms(1,3)=0.5  ! Ti atoms
 O_atoms(1,1)=0.5;  O_atoms(1,2)=0.5;  O_atoms(1,3)=0.0  ! O1 atoms ! xy plane
 O_atoms(2,1)=0.5;  O_atoms(2,2)=0.0;  O_atoms(2,3)=0.5  ! O2 atoms ! xz plane
 O_atoms(3,1)=0.0;  O_atoms(3,2)=0.5;  O_atoms(3,3)=0.5  ! O3 atoms ! yz plane
!!!!!!!!!!!!!!!!!!!!!!!printing input files!!!!!!!!
PRINT*, " *********************************"
PRINT*, "!!!!!!! ABO3 unit cell !!!!!!!"
       DO i=1,3
        PRINT*, a_unit(i,1), a_unit(i,2), a_unit(i,3)
        ENDDO

        PRINT*, "Ba atom coordinate", Ba_atoms(1,1), Ba_atoms(1,2), Ba_atoms(1,3)
        PRINT*, "Ti atom coordinate", Ti_atoms(1,1),Ti_atoms(1,2),Ti_atoms(1,3)
        PRINT*, " O atoms coordinate"
         DO i=1,3
           PRINT*, O_atoms(i,1),O_atoms(i,2),O_atoms(i,3)
           ENDDO
PRINT*, "*************END********************"
!a_unit_T=TRANSPOSE(a_unit)
WRITE(13,*) " ABO3  "
WRITE(13,*) 1.0
WRITE(13,300)  a_strain(1,1)*nx, a_strain(1,2),   a_strain(1,3)
WRITE(13,300)  a_strain(2,1),    a_strain(2,2)*ny,a_strain(2,3)
WRITE(13,300)  a_strain(3,1),    a_strain(3,2),   a_strain(3,3)*nz

WRITE(13,*) "Ba  Ti  O"
WRITE(13,301) ntot, ntot, 3*ntot
301 FORMAT(3I7)
WRITE(13,*) "Cartesian"        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Supercell!!!!!!!!!!!!!!!!!!!!!!!
 Ba_o(1,1)=Ba_atoms(1,1)*a_unit(1,1)
 Ba_o(1,2)=Ba_atoms(1,2)*a_unit(2,2)
 Ba_o(1,3)=Ba_atoms(1,3)*a_unit(3,3)
 Ti_o(1,1)=Ti_atoms(1,1)*a_unit(1,1)
 Ti_o(1,2)=Ti_atoms(1,2)*a_unit(2,2)
 Ti_o(1,3)=Ti_atoms(1,3)*a_unit(3,3)


 Ba_s(1,1)=Ba_atoms(1,1)*a_strain(1,1)
 Ba_s(1,2)=Ba_atoms(1,2)*a_strain(2,2)
 Ba_s(1,3)=Ba_atoms(1,3)*a_strain(3,3)

 Ti_s(1,1)=Ti_atoms(1,1)*a_strain(1,1)
 Ti_s(1,2)=Ti_atoms(1,2)*a_strain(2,2)
 Ti_s(1,3)=Ti_atoms(1,3)*a_strain(3,3)

!PRINT*, Ti_s(1,1), Ti_s(1,2), Ti_s(1,3)

!!!!!Ba atoms
n=0
          DO k=0, nz-1 ! z-direction 
            DO j=0, ny-1 ! y-direction
              DO i=0, nx-1 ! x-direction
               n=n+1
               T=(/i*a_unit(1,1), j*a_unit(2,2), k*a_unit(3,3)/) 
               T_s=(/i*a_strain(1,1), j*a_strain(2,2), k*a_strain(3,3)/)

               b=(/T(1)+Ba_o(1,1), T(2)+Ba_o(1,2), T(3)+Ba_o(1,3)/)
               A_atoms(n,:)=(/b(1), b(2), b(3)/)
               !PRINT*, T_s(1), T_s(2), T_s(3)
               b=(/T_s(1)+Ba_s(1,1), T_s(2)+Ba_s(1,2), T_s(3)+Ba_s(1,3)/)
               CALL domain_struc_A(i, b, bf)
               A_s(n,:)=(/bf(1), bf(2), bf(3)/)
               WRITE(13,300) bf(1), bf(2), bf(3)
                ENDDO
                  ENDDO
                      ENDDO
!!!! Ti atoms
          DO k=0, nz-1
            DO j=0, ny-1
              DO i=0, nx-1
               n=n+1
               T=(/i*a_unit(1,1), j*a_unit(2,2), k*a_unit(3,3)/)
               T_s=(/i*a_strain(1,1), j*a_strain(2,2), k*a_strain(3,3)/)

               b=(/T(1)+Ti_o(1,1), T(2)+Ti_o(1,2), T(3)+Ti_o(1,3)/)
               B_atoms(n,:)=(/b(1), b(2), b(3)/)
               !PRINT*, T_s(1), T_s(2), T_s(3)
            
               b=(/T_s(1)+Ti_s(1,1), T_s(2)+Ti_s(1,2), T_s(3)+Ti_s(1,3)/)

               CALL domain_struc_B(i, b, bf)
               B_s(n,:)=(/bf(1), bf(2), bf(3)/)
                !PRINT*, bf(1), bf(2), bf(3)   
                WRITE(13,300) bf(1), bf(2), bf(3)
                ENDDO
                  ENDDO
                      ENDDO
!!!O1 atoms!!!!
 O_o(1,1)=O_atoms(1,1)*a_unit(1,1)
 O_o(1,2)=O_atoms(1,2)*a_unit(2,2)
 O_o(1,3)=O_atoms(1,3)*a_unit(3,3) 
  
 O_s(1,1)=O_atoms(1,1)*a_strain(1,1)
 O_s(1,2)=O_atoms(1,2)*a_strain(2,2)
 O_s(1,3)=O_atoms(1,3)*a_strain(3,3)
      
          DO k=0, nz-1
            DO j=0, ny-1
              DO i=0, nx-1
               n=n+1
               T=(/i*a_unit(1,1), j*a_unit(2,2), k*a_unit(3,3)/)
               T_s=(/i*a_strain(1,1), j*a_strain(2,2), k*a_strain(3,3)/)

               b=(/T(1)+O_o(1,1), T(2)+O_o(1,2), T(3)+O_o(1,3)/)
               O1_atoms(n,:)=(/b(1), b(2), b(3)/)

               b=(/T_s(1)+O_s(1,1), T_s(2)+O_s(1,2), T_s(3)+O_s(1,3)/)
               CALL domain_struc_O(i, b, bf)
               O1_s(n,:)=(/bf(1), bf(2), bf(3)/)
               WRITE(13,300) bf(1), bf(2), bf(3)
              ! WRITE(13,*) b(1), b(2), b(3)
              ! PRINT*, n, b(1), b(2), b(3)
                ENDDO
                  ENDDO
                      ENDDO
!!!!!!!!O2 atoms
 O_o(2,1)=O_atoms(2,1)*a_unit(1,1)
 O_o(2,2)=O_atoms(2,2)*a_unit(2,2)
 O_o(2,3)=O_atoms(2,3)*a_unit(3,3)

 O_s(2,1)=O_atoms(2,1)*a_strain(1,1)
 O_s(2,2)=O_atoms(2,2)*a_strain(2,2)
 O_s(2,3)=O_atoms(2,3)*a_strain(3,3)

          DO k=0, nz-1
            DO j=0, ny-1
              DO i=0, nx-1
               n=n+1
               T=(/i*a_unit(1,1), j*a_unit(2,2), k*a_unit(3,3)/)
               T_s=(/i*a_strain(1,1), j*a_strain(2,2), k*a_strain(3,3)/)

               b=(/T(1)+O_o(2,1), T(2)+O_o(2,2), T(3)+O_o(2,3)/)
               O2_atoms(n,:)=(/b(1), b(2), b(3)/)
               
              b=(/T_s(1)+O_s(2,1), T_s(2)+O_s(2,2), T_s(3)+O_s(2,3)/)
              CALL domain_struc_O(i, b, bf)
              O2_s(n,:)=(/bf(1), bf(2), bf(3)/)
              WRITE(13,300) bf(1), bf(2), bf(3)
               !WRITE(13,*) b(1), b(2), b(3)
                ENDDO
                  ENDDO
                      ENDDO

!!!!!!O3 atoms
 O_o(3,1)=O_atoms(3,1)*a_unit(1,1)
 O_o(3,2)=O_atoms(3,2)*a_unit(2,2)
 O_o(3,3)=O_atoms(3,3)*a_unit(3,3)

 O_s(3,1)=O_atoms(3,1)*a_strain(1,1)
 O_s(3,2)=O_atoms(3,2)*a_strain(2,2)
 O_s(3,3)=O_atoms(3,3)*a_strain(3,3)

          DO k=0, nz-1
            DO j=0, ny-1
              DO i=0, nx-1
               n=n+1
               T=(/i*a_unit(1,1), j*a_unit(2,2), k*a_unit(3,3)/)
               T_s=(/i*a_strain(1,1), j*a_strain(2,2), k*a_strain(3,3)/)

               b=(/T(1)+O_o(3,1), T(2)+O_o(3,2), T(3)+O_o(3,3)/)
               O3_atoms(n,:)=(/b(1), b(2), b(3)/)
               
               b=(/T_s(1)+O_s(3,1), T_s(2)+O_s(3,2), T_s(3)+O_s(3,3)/)
               CALL domain_struc_O(i, b, bf)
               O3_s(n,:)=(/bf(1), bf(2), bf(3)/)
               WRITE(13,300) bf(1), bf(2), bf(3)
!              PRINT*, n
               !WRITE(13,*) b(1), b(2), b(3)
                ENDDO
                  ENDDO
                      ENDDO
300  FORMAT(3F12.7)
!!!!!!!!!!!!!!!!!!!! Writing lammps dada file format!!!!!!!!!
WRITE(15,*)'#LAMMPS ABO3 domain'
WRITE(15,*)"  "
WRITE(15,190) 5*ntot, "atoms"
WRITE(15,191) 3*ntot, "angles"

PRINT*, " Supercell(nx*ny*nz)=", nx, ny, nz
PRINT*, " Total number of atoms=", 5*ntot
PRINT*,"**************Warning:********************"
PRINT*, "If total number of atoms greater than ten million atoms, Please modify FORMAT and INTEGER options" 
PRINT*, "                               "

190 FORMAT(I9, A6)
191 FORMAT(I9, 1X, A6)

WRITE(15,192) 3, "atom", "types"
WRITE(15,192) 1, "angle", "types"
WRITE(15,*)"  "

192 FORMAT(I2, A6, A6)

WRITE(15,200)0.0, nx*a_strain(1,1),  "xlo xhi"
WRITE(15,200)0.0, ny*a_strain(2,2),  "ylo yhi"
WRITE(15,200)0.0, nz*a_strain(3,3),  "zlo zhi"
200   FORMAT(1X, F3.1, F12.5, 2X, A7)

WRITE(15,*)"  "
WRITE(15,*)"Masses"
WRITE(15,*)"  "

PRINT*,"**************Warning:********************"
PRINT*, "if the structure is not BaTiO3 Please modify atomic masses" 
PRINT*, "                                             "

WRITE(15,202)"1", massA ! " # Ba" !Modify this for other ABO3 compounds
WRITE(15,202)"2", massB !" # Ti"   ! Modify this for other ABO3 compounds
WRITE(15,202)"3", massO !" # O" ! Modify this for other ABO3 compounds
WRITE(15,202)"   "
201 FORMAT(I2, F9.4)
202 FORMAT(A2, F12.5)

WRITE(15,*)"Angle Coeffs"
WRITE(15,*)"  "

PRINT*,"**************Warning:********************"
PRINT*, "if the structure is not BaTiO3 Please modify Angle Coeffs "
PRINT*, "                                                        "

WRITE(15,205)1, angle_coeffs,  180.00 ! Modify this for other ABO3 compounds
205 FORMAT(1X,I1, 1X, I3, 1X, F7.2)

WRITE(15,*)"  "
WRITE(15,*)"Atoms"
WRITE(15,*)"  "
!Modify this for other ABO3 compounds

PRINT*,"**************Warning:********************"
PRINT*, "if the structure is not BaTiO3 Please modify atomic charges"
PRINT*, "                                                        "

!charge_A=1.34730 
!charge_B=1.28905
!charge_O=-0.87878
!
  DO i=1, ntot
   WRITE(15,100)i, 1, 1, charge_A, A_s(i, 1), A_s(i, 2), A_s(i, 3)
   !PRINT*, "ok"
    ENDDO
DO i=nb, 2*ntot
!  PRINT*, "ok"
   WRITE(15,100)i, 1, 2, charge_B, B_s(i, 1), B_s(i, 2), B_s(i, 3)
ENDDO

DO i=no1, 3*ntot
  WRITE(15,100)i, 1, 3, charge_O, O1_s(i, 1), O1_s(i, 2), O1_s(i, 3)
  ENDDO

DO i=no2,4*ntot
  WRITE(15,100)i, 1, 3, charge_O, O2_s(i, 1), O2_s(i, 2), O2_s(i, 3)
  ENDDO

DO i=no3, 5*ntot
  WRITE(15,100)i, 1, 3, charge_O, O3_s(i, 1), O3_s(i, 2), O3_s(i, 3)
  ENDDO

100 FORMAT(I10, 1X, I3, 1X, I3, 1X, F10.5, 1X, 3F12.7)

WRITE(15,*)"  "
WRITE(15,*)"Angles"
WRITE(15,*)"  "

!!! Making a list of 180 angle among the oxigen!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!O1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nx=1,2,3,4.... layer along x-axis
! ny=1,2,3,4....layer along y-axis
! nz=1,2,3,4.... lyer along z-axis
! nx*ny is the number of atoms in xy plane
! DO loop making supercell are
! DO i=0,nz-1
!    DO j=0, ny-1
!      DO k=0, nx-1
!
!  nz
!   |
!   |____ nx
!   / 
!  /
! nx
!  
!
m1=0
n=2*ntot
IF (nz>2) THEN
DO i=1, nz-2 ! z-direction layers
  DO j=1, nxy  ! 1st layer in nx*ny plane
   n=n+1
   m1=m1+1
   n1o=2*ntot+(i-1)*nxy+j; n2o=2*ntot+i*nxy+j; n3o=2*ntot+(i+1)*nxy+j
   CALL ANGLE_AMONG1(O1_atoms(n1o,:),O1_atoms(n2o,:), O1_atoms(n3o,:), angle )
   WRITE(14,*) m1, 1, n1o, n2o, n3o, angle
   WRITE(15,*)m1, 1, n1o, n2o, n3o
   !WRITE(14,*) n, 1, 2*ntot+(i-1)*nxy+j, 2*ntot+i*nxy+j, 2*ntot+(i+1)*nxy+j
   !PRINT*, 2*ntot+(i-1)*nxy+j, 2*ntot+i*nxy+j, 2*ntot+(i+1)*nxy+j
   ENDDO
    ENDDO
ENDIF
!! nxy layer 
  DO j=1, nxy
   n=n+1
   m1=m1+1
    n1o=2*ntot+(nz-1)*nxy+j; n2o=2*ntot+ j; n3o=2*ntot+ nxy+j
    CALL ANGLE_AMONG2(O1_atoms(n1o,:),O1_atoms(n2o,:), O1_atoms(n3o,:), angle )
    WRITE(14,*) m1, 1, n1o, n2o, n3o, angle
    WRITE(15,*)m1, 1, n1o, n2o, n3o
   ENDDO
!nxy-1 layer
 DO j=1, nxy
  n=n+1
  m1=m1+1
  n1o=2*ntot+(nz-2)*nxy+j; n2o=2*ntot+(nz-1)*nxy+j; n3o=2*ntot+j
 ! WRITE(14,*) n, 1,  2*ntot+(nz-2)*nxy+j, 2*ntot+(nz-1)*nxy+j, 2*ntot+j
  CALL ANGLE_AMONG2(O1_atoms(n1o,:),O1_atoms(n2o,:), O1_atoms(n3o,:), angle )
  WRITE(14,*) m1, 1, n1o, n2o, n3o, angle
  WRITE(15,*)m1, 1, n1o, n2o, n3o
  ENDDO

PRINT*, "Number of angle for O1 atoms=", m1

!!!!!!!!!!!!!!!!!!!O2!!!!!!!!!!!!!!!!
n=3*ntot
!IF (nz>2) THEN
DO k=0, nz-1
  DO i=1, ny-2 ! y-direction layer
    DO j=1, nx
     n=n+1
     m1=m1+1
     n1o=3*ntot+(i-1)*nx+j + k*nxy; n2o=3*ntot+i*nx+j+k*nxy; n3o=3*ntot+(i+1)*nx+j+k*nxy
     CALL ANGLE_AMONG1(O2_atoms(n1o,:),O2_atoms(n2o,:), O2_atoms(n3o,:), angle )
     WRITE(14,*) m1, 1, n1o, n2o, n3o, angle
     WRITE(15,*)m1, 1, n1o, n2o, n3o
     !WRITE(14,*) n, 1, 3*ntot+(i-1)*nx+j + k*nxy , 3*ntot+i*nx+j+k*nxy, 3*ntot+(i+1)*nx+j+k*nxy
     ENDDO
      ENDDO
        ENDDO
!ENDIF
!!!!!
!!!!!
DO k=0, nz-1
  DO j=1, nx
   n=n+1
   m1=m1+1
   n1o=3*ntot+nx+j + k*nxy; n2o=3*ntot+j+k*nxy; n3o=3*ntot+(ny-1)*nx+j+k*nxy
   CALL ANGLE_AMONG2(O2_atoms(n1o,:),O2_atoms(n2o,:), O2_atoms(n3o,:), angle)
   WRITE(14,*) m1, 1, n1o, n2o, n3o, angle
   WRITE(15,*)m1, 1, n1o, n2o, n3o
   !WRITE(14,*) n, 1, 3*ntot+nx+j + k*nxy , 3*ntot+j+k*nxy, 3*ntot+(ny-1)*nx+j+k*nxy
   ENDDO
     ENDDO

!
DO k=0, nz-1
  DO j=1, nx
  n=n+1
  m1=m1+1
  n1o=3*ntot+(ny-2)*nx+j + k*nxy; n2o=3*ntot+(ny-1)*nx+j+ k*nxy; n3o=3*ntot+j+ k*nxy
  CALL ANGLE_AMONG2(O2_atoms(n1o,:),O2_atoms(n2o,:), O2_atoms(n3o,:), angle)
  WRITE(14,*) m1, 1, n1o, n2o, n3o, angle
  WRITE(15,*)m1, 1, n1o, n2o, n3o
   !WRITE(14,*) n, 1, 3*ntot+(ny-2)*nx+j + k*nxy , 3*ntot+(ny-1)*nx+j+ k*nxy, 3*ntot+j+ k*nxy
   ENDDO
     ENDDO
PRINT*, "Number of angle for O2 atoms=", m1-ntot

!------------------------------ O3----------------------------
n=4*ntot
!PRINT*, m1

IF (nx>2) THEN
DO i=0, ny-1
DO k=0, nz-1
  DO j=1, nx-2
n=n+1
m1=m1+1
 n1o=4*ntot+j + (i*nx) + k*nxy;n2o=4*ntot+(j+1)+ (i*nx) + k*nxy;n3o=4*ntot+(j+2)+(i*nx)+ nxy*k
  CALL ANGLE_AMONG1(O3_atoms(n1o,:),O3_atoms(n2o,:), O3_atoms(n3o,:), angle)
!  PRINT*, n1o, n2o, n3o
  WRITE(14,*) m1, 1, n1o, n2o, n3o, angle
  WRITE(15,*)m1, 1, n1o, n2o, n3o
  ! WRITE(13,*) n, 1, 4*ntot+j + i*nx + k*nxy , 4*ntot+j+1+ i*nx+ k*nxy, 4*ntot+j+2+i*nx+ nxy*k
   ENDDO
     ENDDO
        ENDDO
ENDIF

!PRINT*, m1

DO k=0, nz-1
  DO i=0, ny-1
n=n+1
m1=m1+1
 n1o=4*ntot+ 2 + (i*nx) + k*nxy ;n2o=4*ntot+ 1 + (i*nx) + k*nxy ;n3o=4*ntot+ (i*nx) + nx + nxy*k 
  CALL ANGLE_AMONG2(O3_atoms(n1o,:),O3_atoms(n2o,:), O3_atoms(n3o,:), angle)
!  PRINT*, n1o, n2o, n3o
  WRITE(14,*) m1, 1, n1o, n2o, n3o, angle
  WRITE(15,*)m1, 1, n1o, n2o, n3o
!
   !WRITE(13,*) n, 1, 4*ntot+2 + i*nx + k*nxy , 4*ntot+ 1 + i*nx+ k*nxy, 4*ntot+i*nx+ 10 + nxy*k
   ENDDO
     ENDDO

!PRINT*, m1


IF(nx<=2) THEN

DO j=1, ny-2
DO k=0, nz-1
  DO i=2, nx
n=n+1
m1=m1+1
 n1o=4*ntot+ (2*nx) + k*nxy + j*nx ; n2o=4*ntot + (2*nx) + k*nxy + (j*nx-1) ;n3o=4*ntot+ (2*nx) +  nxy*k + j*nx
  CALL ANGLE_AMONG2(O3_atoms(n1o,:),O3_atoms(n2o,:), O3_atoms(n3o,:), angle)
 ! PRINT*, n1o, n2o, n3o
  WRITE(14,*) m1, 1, n1o, n2o, n3o, angle
  WRITE(15,*)m1, 1, n1o, n2o, n3o
!
   !WRITE(13,*) n, 1, 4*ntot+2 + i*nx + k*nxy , 4*ntot+ 1 + i*nx+ k*nxy, 4*ntot+i*nx+ 10 + nxy*k
   ENDDO
     ENDDO
ENDDO
ENDIF


DO k=0, nz-1  !nz=ny
  DO i=0, ny-1
n=n+1
m1=m1+1
 n1o=4*ntot+ (nx-1) + (i*nx) + k*nxy; n2o=4*ntot+ nx + (i*nx) + k*nxy; n3o=4*ntot+ (i*nx) + 1 + nxy*k
  CALL ANGLE_AMONG2(O3_atoms(n1o,:),O3_atoms(n2o,:), O3_atoms(n3o,:), angle)
!  PRINT*, n1o, n2o, n3o
  WRITE(14,*) m1, 1, n1o, n2o, n3o, angle
  WRITE(15,*)m1, 1, n1o, n2o, n3o
!   !WRITE(13,*) n, 1, 4*ntot+9 + i*nx + k*nxy , 4*ntot+ 10 + i*nx+ k*nxy, 4*ntot+i*nx+ 1 + nxy*k
   ENDDO
     ENDDO

IF (nx<=2) THEN
DO j=1, ny-2
DO k=0, nz-1
  DO i=2, nx
n=n+1
m1=m1+1

 n1o=4*ntot+ (j*nx-1) + k*nxy + 2*nx; n2o=4*ntot+ (2*nx) + k*nxy +j*nx ; n3o=4*ntot + (j*nx-1) + 2*nx + nxy*k
  CALL ANGLE_AMONG2(O3_atoms(n1o,:),O3_atoms(n2o,:), O3_atoms(n3o,:), angle)
  !PRINT*, n1o, n2o, n3o
  WRITE(14,*) m1, 1, n1o, n2o, n3o, angle
  WRITE(15,*)m1, 1, n1o, n2o, n3o
!
   !WRITE(13,*) n, 1, 4*ntot+2 + i*nx + k*nxy , 4*ntot+ 1 + i*nx+ k*nxy, 4*ntot+i*nx+ 10 + nxy*k
   ENDDO
     ENDDO
ENDDO
ENDIF
PRINT*, "Number of angle for O3 atoms=", m1-2*ntot
PRINT*, "                                      "
PRINT*, "************ Warning:*****************"
PRINT*, " Please check the output of this program before the final calculation"
PRINT*, " Their is no grante that the program ouput is always right"
PRINT*,"                                       "   
!!
         END PROGRAM ABO3_domain_struc
