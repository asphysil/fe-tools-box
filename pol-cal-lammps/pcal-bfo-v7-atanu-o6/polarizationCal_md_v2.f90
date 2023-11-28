! Purpose to calculate polarization from a structure

PROGRAM POLARIZATION_CAL
use user_fileID 
use constant
use data_structure
use nearest_neighbour_list, only : nearest_neighbour_ab, nearest_neighbour_aa,&
nearest_neighbour_ab1, nearest_neighbour_aa1

use relative_to_O_disp, only :  cal_a_site_pol_disp,cal_b_site_pol_disp

use writefiles, only : write_xsf, write_pol,&
               write_dump_ovito, write_pol_ovito, &
               write_dump_ovito_LMP

use readfiles, only : read_dumpxyz_spin,&
                      read_angle_info_spin
IMPLICIT NONE   
!     REAL(dp), DIMENSION(:,:) :: frac_pb(3,na),frac_ti(3,na), frac_O(3,3*na)

     REAL(dp),DIMENSION(ndim3) ::  ptemp

     REAL(dp) :: x,y, start, finish, &
                vol,  nti_dp, born_charge
    REAL(dp) :: dab, daa, dbb  

     INTEGER :: ntot,  npb, nti, no, ndump, tot_ndump,  no_neigh
     INTEGER :: i, j, k, m, junk1_int

     INTEGER :: fileID
     CHARACTER(LEN=40):: filename
     LOGICAl :: fexist

! Initialization variables
npb=0
nti=0
no=0
DO j=1,3
   DO i=1,3
     latt_vec(i,j)=0.0
ENDDO
ENDDO
! File unit No.
DO i = 1, nmax_file
  fileID_all(i) = 40 + i 
ENDDO

!
! reading few lines of dump.xyz file
!PRINT*, " Enter total number of times dumped the atomic coordinates"
!READ*, tot_ndump
tot_ndump = 1
PRINT*,"***************Warning*****************************"
PRINT*,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
PRINT*,"%%%%%%  %%%%%%"
PRINT*, " Please check input values of A-B, A-O, B-O, B-B distance are correctly set"
PRINT*," "
PRINT*," "
PRINT*,"! "
PRINT*,"*********************************************"

PRINT*, " Enter total number of times dumped the atomic coordinates"
READ*, tot_ndump
!
!PRINT*, "Please change total number of times dumped . default is one"
!tot_ndump=1
!!!!!!!!!!!!!!!!!
INQUIRE(FILE="dump.xyz", exist=fexist)
IF (fexist) THEN
    PRINT*, "dump.xyz file exist"
ELSE
  PRINT*, "dump.xyz file does exist in this folder and program will STOP here"
STOP
ENDIF

INQUIRE(FILE="data-jiahao.BFO", exist=fexist)
IF (fexist) THEN
    PRINT*, "data-jiahao.BFO file exist"
ELSE
  PRINT*, "data-jiahao.BFO file does exist in this folder and program will STOP here"
STOP
ENDIF


OPEN(UNIT=12, FILE='dump.xyz', ACTION='READ')
READ(12,*)
DO i=1, 2
  READ(12,*)
ENDDO
READ(12,*) ntot
close(12)
!!!!!!!!!!!!!!!
npb=ntot/5
nti=npb 
no=3*nti
nti_dp =  DBLE(nti)
!==================
! ALLOCATE array 
! Atoms
ALLOCATE(a_atms(ndim3, npb),b_atms(ndim3, npb))
ALLOCATE(o_atms(ndim3, no))
! DISP.
ALLOCATE( a_disp(ndim3, npb), b_disp(ndim3, npb))
! POL.
ALLOCATE( a_pol(ndim3, npb), b_pol(ndim3, npb), pol_unit(ndim3, npb)) 
ALLOCATE(pol_ab(ndim3,no_ab_neigh,npb))
!
!ALLOCATE( a_pol_dir(npb))

! SPIN
ALLOCATE(sp(ndim3, npb))
! Nearest neighbour and angle
ALLOCATE(coord_fe1fe6(ndim3, no_aa_neigh, npb))
ALLOCATE(coord_fe1ofe2(ndim3, ndim3, no))
ALLOCATE(angle_index(ndim3, no))


! Nearest neighbour list
ALLOCATE(aa_neigh(no_aa_neigh, npb))
ALLOCATE(bb_neigh(no_bb_neigh, npb))

ALLOCATE(ab_neigh(no_ab_neigh, npb))
ALLOCATE(ao12_neigh(no_ao_neigh, npb))
ALLOCATE(bo6_neigh(no_bo_neigh, npb))



! openn file 1D 
fileID_read = fileID_all(1)
fileID_totpol = fileID_all(2)
fileID_locpol = fileID_all(3)
fileID_ovito = fileID_all(4)
fileID_pol_ovito = fileID_all(5)
!fileID_rtheta = fileID_all(6)
file_ID_xsf1 = fileID_all(7)
file_ID_xsf2 = fileID_all(8)

fileID_angle = fileID_all(9)

fileID_ovito_lmp =  fileID_all(10)

OPEN(UNIT=fileID_read, FILE='dump.xyz', ACTION='READ')
OPEN(UNIT=fileID_angle, FILE="data-jiahao.BFO", ACTION='READ')

OPEN(UNIT=fileID_totpol,FILE="Total_P.dat", ACTION="WRITE")   
OPEN(UNIT=fileID_locpol,FILE="Local_P.dat", ACTION="WRITE")
OPEN(UNIT=fileID_ovito,FILE="dump-ovito.xyz", ACTION="WRITE")
OPEN(UNIT=fileID_pol_ovito, FILE="dump-pol-ovito.xyz", ACTION="WRITE")
!OPEN(UNIT=fileID_rtheta, FILE="rtheta-plt.dat", ACTION="WRITE")
OPEN(UNIT=fileID_ovito_lmp, FILE="dump-ovito-LMP.xyz", ACTION="WRITE")


! reading dump file
CALL CPU_TIME(start)

CALL read_dumpxyz_spin(fileID_read,npb, nti, no)
CALL read_angle_info_spin(fileID_angle,npb, no)
close(fileID_angle)


CALL CPU_TIME(finish)
PRINT*, "  //////////////// "
PRINT*, " Total CPU time to read a single dump =", (finish-start)/60.0, "Minutes"
         !  cell volume
vol = ABS( latt_vec(1,1) * (latt_vec(2,2)*latt_vec(3,3) - latt_vec(2,3)*latt_vec(3,2)) &
          -latt_vec(2,1) * (latt_vec(1,2)*latt_vec(3,3) - latt_vec(1,3)*latt_vec(3,2)) &
          +latt_vec(3,1) * (latt_vec(1,2)*latt_vec(2,3) - latt_vec(1,3)*latt_vec(2,2)) )        
unitvol=vol/nti
!!for calculating inverse Matrix
DO j=1, 3
    DO i=1, 3
     ainv(i,j)=latt_vec(i,j)
   ENDDO
 ENDDO


!! call DGETRF(ndim, ndim, ainv, ndim, ipiv, info)
!! IF (info /= 0) THEN
!!    PRINT*, 'Matrix is numerically singular!'
!!    STOP
!! ENDIF
!!
!!   call DGETRI(ndim, ainv, ndim, ipiv, work, ndim, info)
!!
!!     IF (info /= 0) THEN
!!        PRINT*, 'Matrix inversion failed!'
!!        STOP
!!     ENDIF
!!
!!! Cartesian to fractional coordinate transformation
!     frac_ti=MATMUL(tib(1:nti,:),ainv)
!     frac_pb=MATMUL(pba(1:npb,:),ainv)
!     frac_O=MATMUL(O(1:no,:),ainv)

!!!! finding nearest neighbour distance
CALL CPU_TIME(start)
INQUIRE(FILE="aa-neighbours.dat", exist=fexist)
IF (fexist) THEN
    PRINT*, "aa-neighbours.dat file exist"

    OPEN(UNIT=12, FILE="aa-neighbours.dat", ACTION="READ")
     do i=1, npb
     READ(12, *) junk1_int, (aa_neigh(j,i), j=1,6)
    enddo
ELSE
  PRINT*, "aa-neighbours.dat file does not exist"
  OPEN(UNIT=12, FILE="aa-neighbours.dat", ACTION="WRITE")
  ! robust subroutine but slow
 ! daa = daa_cut  
 ! no_neigh=no_aa_neigh
 ! CALL  nearest_neighbour_aa(no_neigh, npb, npb, daa, a_atms(:,1:npb), a_atms(:,1:npb), aa_neigh(:,1:npb))
!!
  daa = daa_min 
  no_neigh=no_aa_neigh
  CALL  nearest_neighbour_aa1(no_neigh, npb, npb, daa, a_atms(:,1:npb), a_atms(:,1:npb), aa_neigh(:,1:npb))
   do i=1, npb
   WRITE(12, *) i, (aa_neigh(j,i), j=1,6)
   enddo
 ENDIF
close(12)

!========bb=================
INQUIRE(FILE="bb-neighbours.dat", exist=fexist)
IF (fexist) THEN
    PRINT*, "bb-neighbours.dat file exist"

    OPEN(UNIT=12, FILE="bb-neighbours.dat", ACTION="READ")
     do i=1, nti
     READ(12, *) junk1_int, (bb_neigh(j,i), j=1,6)
    enddo
ELSE
  PRINT*, "bb-neighbours.dat file does not exist"
  OPEN(UNIT=12, FILE="bb-neighbours.dat", ACTION="WRITE")
  ! robust subroutine but slow
 ! daa = daa_cut  
 ! no_neigh=no_aa_neigh
 ! CALL  nearest_neighbour_aa(no_neigh, npb, npb, daa, a_atms(:,1:npb), a_atms(:,1:npb), aa_neigh(:,1:npb))
!!
  dbb = dbb_min 
  no_neigh=no_bb_neigh
  CALL  nearest_neighbour_aa1(no_neigh, nti, nti, daa, b_atms(:,1:nti), b_atms(:,1:nti), bb_neigh(:,1:nti))
   do i=1, nti
   WRITE(12, *) i, (bb_neigh(j,i), j=1,6)
   enddo
 ENDIF
close(12)


!=================


INQUIRE(FILE="ab-neighbours.dat", exist=fexist)
IF (fexist) THEN
    PRINT*, "ab-neighbours.dat file exist"
    
    OPEN(UNIT=13, FILE="ab-neighbours.dat", ACTION="READ")
    do i=1, nti
     READ(13, *) junk1_int, (ab_neigh(j,i), j=1,8)
    enddo
ELSE

  PRINT*, "ab-neighbours.dat file does not exist"
  OPEN(UNIT=13, FILE="ab-neighbours.dat", ACTION="WRITE")
  

  !dab = dab_cut 
  !no_neigh=no_ab_neigh
  !CALL  nearest_neighbour_ab(no_neigh, nti, npb, dab, b_atms(:,1:nti), a_atms(:,1:npb), ab_neigh(:,1:nti))


  dab = dab_min 
  no_neigh=no_ab_neigh
  CALL  nearest_neighbour_ab1(no_neigh, nti, npb, dab, b_atms(:,1:nti), a_atms(:,1:npb), ab_neigh(:,1:nti))
  
   do i=1, nti
   WRITE(13, *) i, (ab_neigh(j,i), j=1,8)
   enddo
  
 ENDIF
close(13)

INQUIRE(FILE="ao12-neighbours.dat", exist=fexist)
IF (fexist) THEN
   PRINT*, "ao12-neighbours.dat file exist"
   
   OPEN(UNIT=14, FILE="ao12-neighbours.dat", ACTION="READ")
   
   do i=1, npb
    READ(14, *) junk1_int, (ao12_neigh(j,i), j=1,12)
   enddo

ELSE
PRINT*, "ao12-neighbours.dat file does not exist"
OPEN(UNIT=14, FILE="ao12-neighbours.dat", ACTION="WRITE")

!dab = dao_cut
!no_neigh=no_ao_neigh
!CALL  nearest_neighbour_ab(no_neigh, npb, no, dab, a_atms(:,1:npb), o_atms(:,1:no),ao12_neigh(:,1:npb))
 
dab = dao_min
no_neigh=no_ao_neigh
CALL  nearest_neighbour_ab1(no_neigh, npb, no, dab, a_atms(:,1:npb), o_atms(:,1:no),ao12_neigh(:,1:npb))
 
do i=1, npb
 WRITE(14, *) i, (ao12_neigh(j,i), j=1,12)
 enddo

 ENDIF

close(14)

INQUIRE(FILE="bo6-neighbours.dat", exist=fexist)
IF (fexist) THEN
   PRINT*, "bo6-neighbours.dat file exist"

   OPEN(UNIT=15, FILE="bo6-neighbours.dat", ACTION="READ")

   do i=1, nti
    READ(15, *) junk1_int, (bo6_neigh(j,i), j=1,6)
   enddo

ELSE
PRINT*, "bo6-neighbours.dat file does not exist"
OPEN(UNIT=15, FILE="bo6-neighbours.dat", ACTION="WRITE")


!dab = dbo_cut
!no_neigh=no_bo_neigh
!CALL  nearest_neighbour_ab(no_neigh, nti, no, dab, b_atms(:,1:nti), o_atms(:,1:no), bo6_neigh(:,1:nti)) 

dab = dbo_min
no_neigh=no_bo_neigh
CALL  nearest_neighbour_ab1(no_neigh, nti, no, dab, b_atms(:,1:nti), o_atms(:,1:no), bo6_neigh(:,1:nti)) 
 do i=1, nti
 WRITE(15, *) i, (bo6_neigh(j,i), j=1,6)
 enddo
 ENDIF
close(15)

!print*, 'ok'
CALL CPU_TIME(finish)
PRINT*, "  //////////////// "
PRINT*, " Total CPU time for calculating nearest neighbour and distance =", (finish-start)/60.0, "Minutes"
!
! Calculate a-site polarization
born_charge=zqa
CALL cal_a_site_pol_disp(born_charge, npb, no_ao_neigh)



! Calculate b-site polarization  
born_charge=zqb
CALL cal_b_site_pol_disp(born_charge, nti,  no_bo_neigh)
 
! arranging data
DO i=1, nti
  DO k=1, no_ab_neigh
    m = ab_neigh(k, i)
   pol_ab(1:3,k,i)= a_pol(1:3,m)
  ENDDO 
ENDDO

  !!Polarization per unit cell
  DO i=1, nti
    ptemp=(/0.0,0.0,0.0/)
    DO k=1, no_ab_neigh
      DO j= 1, 3 
        ptemp(j) = ptemp(j) +  pol_ab(j,k,i)
      ENDDO
    ENDDO ! x y z
    pol_unit(1:3,i) = (b_pol(1:3,i) + (ptemp(1:3)/8.0))
    !print*, (pol_unit(j,i), j=1,3)
ENDDO

!
!CALL cal_afe_order(npb, no_aa_neigh)

! Write  xsf file

filename='xcrysden-init.xsf'
CALL write_xsf(file_ID_xsf1, filename,  npb, ntot)
close(file_ID_xsf1)

! calling to write Total P and local P in a file
CALL write_pol(fileID_totpol, fileID_locpol, npb) 
!CALL write_dump_ovito_LMP(fileID_ovito_lmp, npb)

!CALL write_pol_ovito_pmwo(fileID_pol_ovito,npb)
CALL write_dump_ovito(fileID_ovito, npb)
CALL write_pol_ovito(fileID_pol_ovito,npb) 


PRINT*, " reading dump >1 "

PRINT*, " $$$$$$$$$ "
PRINT*, " Unit of the Local and Total polarization is in C/m^2"
!PRINT*,  (ptot(j), j = 1, 3), SQRT(ptot(1)**2 + ptot(2)**2 + ptot(3)**2)
CALL CPU_TIME(start)
!!!!!!!! Second itrations!!!!!!!!!!!

DO ndump=2, tot_ndump
! reading dump file
  CALL read_dumpxyz_spin(fileID_read,npb, nti, no)
   !!*******************************!!
            !  cell volume
    vol = ABS( latt_vec(1,1) * (latt_vec(2,2)*latt_vec(3,3) - latt_vec(2,3)*latt_vec(3,2)) &
               - latt_vec(2,1) * (latt_vec(1,2)*latt_vec(3,3) - latt_vec(1,3)*latt_vec(3,2)) &
               + latt_vec(3,1) * (latt_vec(1,2)*latt_vec(2,3) - latt_vec(1,3)*latt_vec(2,2)) )
    unitvol=vol/nti
   ! Calculate a-site polarization
   born_charge=zqa
   CALL cal_a_site_pol_disp(born_charge, npb, no_ao_neigh)
            
  ! Calculate b-site polarization
   born_charge=zqb
   CALL cal_b_site_pol_disp(born_charge, nti, no_bo_neigh)
    
   ! arranging data 
   DO i=1, nti
       DO k=1, no_ab_neigh
         m = ab_neigh(k, i)
       pol_ab(1:3,k,i)= a_pol(1:3,m)
       ENDDO 
   ENDDO
   
   !!Polarization per unit cell
    DO i=1, nti
        ptemp=(/0.0,0.0,0.0/)
        DO k=1, no_ab_neigh
           DO j= 1, 3 
           ptemp(j) = ptemp(j) +  pol_ab(j,k,i)
          ENDDO
        ENDDO ! x y z
    pol_unit(1:3,i) = (b_pol(1:3,i) + (ptemp(1:3)/8.0))
   ENDDO
  
!print*, ' ndump'
CALL write_dump_ovito(fileID_ovito, npb) 
CALL write_pol_ovito(fileID_pol_ovito, npb) 
CALL write_pol(fileID_totpol, fileID_locpol,npb)
!CALL write_dump_ovito_LMP(fileID_ovito_lmp, npb)



    IF ( MOD(ndump,100)==0) THEN
           PRINT*, " number of dumps", ndump, "// in step of 100"
    ENDIF
ENDDO ! ndump

filename='xcrysden-final.xsf'
CALL write_xsf(file_ID_xsf2, filename,  npb, ntot)
close(file_ID_xsf2)

CALL CPU_TIME(finish)
PRINT*, "  //////////////// "
PRINT*, " Total CPU time =", (finish-start)/3600.0, "hrs"


close(fileID_read)
close(fileID_ovito)
close(fileID_pol_ovito)
close(fileID_locpol)
close(fileID_totpol)
close(fileID_ovito_lmp)

END PROGRAM POLARIZATION_CAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
