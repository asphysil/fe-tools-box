!!!!!!!!!!!!!!!!!Main Program!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PROGRAM ABO3_domain_struc

       use, intrinsic :: iso_fortran_env
      use mpi_f08

      USE my_constants, ONLY : dp_real, dp_int, nx, ny, nz
      USE data_structure
      USE data_cubic40,  ONLY : supercell_cubic
     
      IMPLICIT NONE
      INTEGER, PARAMETER :: nunit_max=40
      REAL(dp_real), dimension(ndim3, nunit_max) :: atmu,  atmu_opt

      REAL(dp_real) :: vec_unit(3) 

      INTEGER :: nu, nuba, nuti, nuo1, nuo2, nuo3, nuo 
      !
      INTEGER (dp_int):: nscell, nsba, nsti, nso1, nso2, nso3, nso, ntotal_temp
      INTEGER(dp_int) :: nxy, nxz, nyz
      INTEGER :: flag

      INTEGER(dp_int) :: i, j, n ! integer for loop 
     
   integer(kind=int32) :: ierror ! To control errors in MPI calls
   integer(kind=int32) :: rank   ! Unique number received by each process
   integer(kind=int32) :: num_proc ! Total number of processes


   ! Initialize MPI. This must be the first MPI call
   call MPI_Init(ierror)

   ! Get the number of processes
   call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierror)

   ! Get the individual process rank
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)


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
! 
!! for cubic structure
nu=nunit_max
!CALL crys_struc_cubic(nu, nuba, nuti, nuo1, a_unit, atmu)
CALL crys_struc_cubic_r3c(nu, nuba, nuti, nuo1, a_unit, atmu)
CALL crys_struc_cubic_r3c_opt(nu, nuba, nuti, nuo1, a_unit, atmu_opt)
!


!! for orthorombic structure cut at 110
!CALL crys_struc_ortho(nu, nuba, nuti, nuo1, a_unit, atmu)

! Unit cell
nuo2 = nuo1
nuo3 = nuo1
nuo = nuo1 + nuo2 + nuo3 

!print*, nuti, nuo1, nuo2

! Suppercell
nscell = nx*ny*nz
! Supercell atoms
nsba = nscell*nuba 
nsti = nscell*nuti
nso1 = nscell*nuo1
nso2 = nscell*nuo2
nso3 = nscell*nuo3

nxy = nx*ny*nuti ! xy plane
nxz = nx*nz*nuti ! xz plane
nyz = ny*nz*nuti ! yz plane

nso = nso1 + nso2 + nso3  ! No. of O atoms in suppercell

ALLOCATE(a_atms(4,nsba), b_atms(4, nsti), &
         o1_atms(4,nso1), o2_atms(4,nso2), o3_atms(4,nso3),&
         angle_ti_o_ti(3,nsti), all_angle_ti_o_ti(3, 3*nsti), fe_bfo(nsti))
ALLOCATE(a_atms_opt(4,nsba), b_atms_opt(4,nsti), &
         o1_atms_opt(4,nso1), o2_atms_opt(4,nso2), o3_atms_opt(4, nso3))

DO i=1,3
  Do j=1,3
  latt_sup(i,j) = 0.0
  ENDDO
    ENDDO

!      OPEN(UNIT=12, FILE='POSCAR', ACTION='READ')
      OPEN(UNIT=13, FILE='SUPERCELL.vasp', ACTION='WRITE')
 !     OPEN(UNIT=14, FILE='SUPERCELL_angles.dat', ACTION='WRITE')
 !     OPEN(UNIT=15, FILE='data.BTO', ACTION='WRITE')



!a_strain(1,1) = a_unit(1,1)+ (strainx*a_unit(1,1)*0.01)
!a_strain(2,2) = a_unit(2,2)+ (strainy*a_unit(2,2)*0.01)
!a_strain(3,3) = a_unit(3,3)+ (strainz*a_unit(3,3)*0.01)

!!!!!!!!!!!!!!!!!!!!!!!printing input files!!!!!!!!
PRINT*, " *********************************"
PRINT*, "!!!!!!! ABO3 unit cell !!!!!!!"

       DO i=1,3
        PRINT*, a_unit(i,1), a_unit(i,2), a_unit(i,3)
        ENDDO

        n = 0
        DO i =1, nuba
        n = n + 1
        PRINT*, "Ba atom coordinate", atmu(1,n), atmu(2,n), atmu(3,n)
        ENDDO 

        DO i =1, nuti 
         n = n + 1
        PRINT*, "Ti atom coordinate", atmu(1,n), atmu(2,n),atmu(3,n)
        ENDDO 

        PRINT*, " O atoms coordinate"

        DO i=1, nuo1 
        n = n + 1
           PRINT*, atmu(1,n), atmu(2,n), atmu(3,n)
        ENDDO
    
        DO i=1, nuo2 
            n = n + 1
               PRINT*, atmu(1,n), atmu(2,n), atmu(3,n)
         ENDDO
        DO i=1, nuo3 
            n = n + 1
            PRINT*, atmu(1,n), atmu(2,n), atmu(3,n)
        ENDDO

PRINT*, "*************END********************"
! supercell
latt_sup(1,1) = a_unit(1,1)*DBLE(nx)
latt_sup(2,2) = a_unit(2,2)*DBLE(ny)
latt_sup(3,3) = a_unit(3,3)*DBLE(nz)

! Strain Suppercell
!a_strain(1,1) = a_strain(1,1)*DBLE(nx) 
!a_strain(2,2) = a_strain(2,2)*DBLE(ny)
!a_strain(3,3) = a_strain(3,3)*DBLE(nz)
!
!a_unit_T=TRANSPOSE(a_unit)
WRITE(13,*) " ABO3  "
WRITE(13,*) 1.0

DO i =1, 3
WRITE(13,300)  latt_sup(i,1), latt_sup(i,2),   latt_sup(i,3)
ENDDO 

WRITE(13,*) "Ba  Ti  O"
WRITE(13,301) nsba, nsti, nso 
301 FORMAT(3I7)
WRITE(13,*) "Cartesian"        


flag=1
CALL supercell_cubic(nu, nuti, atmu, flag, ntotal_temp)
flag=2
CALL supercell_cubic(nu, nuti, atmu_opt,flag, ntotal_temp)

DO i = 1, nsba
WRITE(13,*) (a_atms(j, i), j=2,4)
ENDDO

DO i = 1, nsti
WRITE(13,*) (b_atms(j, i), j=2,4)
ENDDO

DO i = 1, nso1
WRITE(13,*) (o1_atms(j, i), j=2,4)
ENDDO
DO i = 1, nso2
WRITE(13,*) (o2_atms(j, i), j=2,4)
ENDDO
DO i = 1, nso3
WRITE(13,*) (o3_atms(j, i), j=2,4)
ENDDO


! Angle along z
vec_unit =(/0.0, 0.0, 1.0/)
CALL interaction_3atms(vec_unit, latt_sup, nsti, b_atms(:,1:nsti), o1_atms(:, 1:nso1), angle_ti_o_ti(:,1:nsti))

  n=0
 DO i = 1, nsti 
    n = n + 1
    all_angle_ti_o_ti(:, n) =  angle_ti_o_ti(:,i)
   !print*, (all_angle_ti_o_ti(n,j), j=1,3)
 ENDDO
 !

!!
! Angle along y
vec_unit =(/0.0, 1.0, 0.0/)
CALL interaction_3atms(vec_unit, latt_sup, nsti,b_atms(:,1:nsti), o2_atms(:, 1:nso2), angle_ti_o_ti(:,1:nsti))
!!PRINT*, '----main-----'
DO i =1,  nsti 
    n = n + 1
    all_angle_ti_o_ti(:, n) =  angle_ti_o_ti(:, i)
   
  !print*, (all_angle_ti_o_ti(n,j), j=1,3)
  !print*, (b_atms(all_angle_ti_o_ti(n,1), j), j=2,4), (o2_atms(all_angle_ti_o_ti(n,2), j), j=2,4)
ENDDO
!

!! ! Angle along x
vec_unit =(/1.0, 0.0, 0.0/)
CALL interaction_3atms(vec_unit, latt_sup, nsti, b_atms(:,1:nsti), o3_atms(:, 1:nso3), angle_ti_o_ti(:,1:nsti))

! !
! 
 DO i = 1, nsti 
     n = n + 1
     all_angle_ti_o_ti(:, n) =  angle_ti_o_ti(:, i)
    !print*, (all_angle_ti_o_ti(n,j), j=1,3)
 ENDDO
!
!
!PRINT*, '----main-----'
!DO i=1, 3*nsti
!    PRINT*, i, ' --- ', (all_angle_ti_o_ti(i,j), j =1,3)
! ENDDO
!

 CALL fe_atoms_afe(nsti, b_atms(:,1:nsti),  fe_bfo(1:nsti) )
! CALL lammps_struc_inp_dmg(nsba, nsti, nso1, nso2, nso3, latt_sup, a_atms, b_atms,&
!                       o1_atms, o2_atms, o3_atms, all_angle_ti_o_ti, fe_bfo)
!DO i=1,nsti
!print*, fe_bfo(i)
!ENDDO

flag=2
CALL  lammps_struc_inp_dmg(flag, nsba, nsti, nso1, nso2, nso3, all_angle_ti_o_ti(:, 1:3*nsti), fe_bfo(1:nsti))

!CALL lammps_struc_inp_liu(nsba, nsti, nso1, nso2, nso3)

!CALL struc_vasp_format(nsba, nsti, nso1, nso2, nso3)


300 FORMAT(3F12.8)

  ! No more MPI calls after Finalize
   call MPI_Finalize(ierror)

END PROGRAM ABO3_domain_struc
