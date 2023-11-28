! Purpose to calculate polarization from a structure

PROGRAM POLARIZATION_CAL
use constant
use data_structure
use nearest_neighbour_list, only : nearest_neighbour_ab, nearest_neighbour_aa,&
nearest_neighbour_ab1, nearest_neighbour_aa1

use relative_to_O_disp, only :  cal_a_site_pol_disp,cal_b_site_pol_disp,&
                                cal_afe_order
use writefiles, only : write_xsf, write_pol, write_dump_ovito, write_pol_ovito,&
                      write_pol_ovito_pmwo, write_pol_afe
use readfiles, only : read_dumpxyz, read_dumpxyz_drew,& 
          readdump_type1, readdump_type2
IMPLICIT NONE   
!     REAL(dp), DIMENSION(:,:) :: frac_pb(3,na),frac_ti(3,na), frac_O(3,3*na)

     REAL(dp),DIMENSION(ndim3) ::  ptemp

     REAL(dp) :: x,y, start, finish, &
                vol,  nti_dp, born_charge
    REAL(dp) :: dab, daa  
    REAL(dp), dimension(ndim3) :: pup, pdown 
    real(dp) :: angle1, norm, inorm  

     INTEGER :: ntot,  npb, nti, no, ndump, tot_ndump,  no_neigh
     INTEGER :: i, j, k, m, junk1_int, m1, m2

     INTEGER :: fileID, fileID_dump
     CHARACTER(LEN=40):: filename
     LOGICAl :: fexist

npb=0
nti=0
no=0


DO j=1,3
   DO i=1,3
     latt_vec(i,j)=0.0
ENDDO
ENDDO

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

! ALLOCATE 
ALLOCATE(a_atms(ndim3, npb),b_atms(ndim3, npb))
allocate(a_atms_symm(ndim3, npb), b_atms_symm(ndim3, npb))

ALLOCATE( a_disp(ndim3, npb), b_disp(ndim3, npb))
ALLOCATE( a_pol(ndim3, npb), b_pol(ndim3, npb), pol_unit(ndim3, npb)) 

ALLOCATE( a_pol_dir(npb))

ALLOCATE(sp(ndim3, npb))

ALLOCATE(o_atms(ndim3, no))
ALLOCATE(pol_ab(ndim3,no_ab_neigh,npb))

ALLOCATE(aa_neigh(no_aa_neigh, npb))

ALLOCATE(ab_neigh(no_ab_neigh, npb))
ALLOCATE(ao12_neigh(no_ao_neigh, npb))
ALLOCATE(bo6_neigh(no_bo_neigh, npb))

allocate(a_pol_up(ndim3, npb), a_pol_down(ndim3, npb))

allocate(pol_up_index(npb/2), pol_down_index(npb/2))

OPEN(UNIT=fileID_in, FILE='dump.xyz', ACTION='READ')

! 
OPEN(UNIT=fileID_totpol,FILE="Total_P.dat", ACTION="WRITE")   
OPEN(UNIT=fileID_locpol,FILE="Local_P.dat", ACTION="WRITE")

OPEN(UNIT=fileID_ovito,FILE="dump-ovito.xyz", ACTION="WRITE")
OPEN(UNIT=fileID_pol_ovito, FILE="dump-pol-ovito.xyz", ACTION="WRITE")
OPEN(UNIT=fileID_rtheta_plt, FILE="rtheta-plt.dat", ACTION="WRITE")
OPEN(UNIT=fileID_afe, FILE="afe-sublattice.dat", ACTION="WRITE")

OPEN(UNIT=fileID_disp_a, FILE="disp-pb-atms.dat", ACTION="WRITE")
OPEN(UNIT=fileID_disp_a_sublatt, FILE="disp-pb-atms-sublatt.dat", ACTION="WRITE")

OPEN(UNIT=fileID_disp_b, FILE="disp-b-site-atom.dat", ACTION="WRITE")
! reading dump file
CALL CPU_TIME(start)
!CALL read_dumpxyz(npb, nti, no)
CALL read_dumpxyz_drew(npb, nti, no)

!CALL readdump_type1(npb, nti, no)
!CALL readdump_type2(npb, nti, no)
!
!fileID=29
!filename='xcrysden-init.xsf'
!CALL write_xsf(fileID, filename,  npb, ntot)

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

    OPEN(UNIT=fileID_in2, FILE="aa-neighbours.dat", ACTION="READ")
    do i=1, npb
     READ(fileID_in2, *) junk1_int, (aa_neigh(j,i), j=1,6)
    enddo
ELSE

  PRINT*, "aa-neighbours.dat file does not exist"
  OPEN(UNIT=fileID_in2, FILE="aa-neighbours.dat", ACTION="WRITE")
  
  ! robust subroutine but slow
  daa = daa_cut  
  no_neigh=no_aa_neigh
  CALL  nearest_neighbour_aa(no_neigh, npb, npb, daa, a_atms(:,1:npb), a_atms(:,1:npb), aa_neigh(:,1:npb))
!!
  
 ! daa = daa_min 
 ! no_neigh=no_aa_neigh
 ! CALL  nearest_neighbour_aa1(no_neigh, npb, npb, daa, a_atms(:,1:npb), a_atms(:,1:npb), aa_neigh(:,1:npb))


   do i=1, npb
   WRITE(fileID_in2, *) i, (aa_neigh(j,i), j=1,6)
   enddo

 ENDIF





INQUIRE(FILE="ab-neighbours.dat", exist=fexist)
IF (fexist) THEN
    PRINT*, "ab-neighbours.dat file exist"
    
    OPEN(UNIT=fileID_in2, FILE="ab-neighbours.dat", ACTION="READ")
    do i=1, nti
     READ(fileID_in2, *) junk1_int, (ab_neigh(j,i), j=1,8)
    enddo
ELSE

  PRINT*, "ab-neighbours.dat file does not exist"
  OPEN(UNIT=fileID_in2, FILE="ab-neighbours.dat", ACTION="WRITE")
  

  dab = dab_cut 
  no_neigh=no_ab_neigh
  CALL  nearest_neighbour_ab(no_neigh, nti, npb, dab, b_atms(:,1:nti), a_atms(:,1:npb), ab_neigh(:,1:nti))


 ! dab = dab_min 
 ! no_neigh=no_ab_neigh
 ! CALL  nearest_neighbour_ab1(no_neigh, nti, npb, dab, b_atms(:,1:nti), a_atms(:,1:npb), ab_neigh(:,1:nti))
  
   do i=1, nti
   WRITE(fileID_in2, *) i, (ab_neigh(j,i), j=1,8)
   enddo
  
 ENDIF


INQUIRE(FILE="ao12-neighbours.dat", exist=fexist)
IF (fexist) THEN
   PRINT*, "ao12-neighbours.dat file exist"
   
   OPEN(UNIT=fileID_in3, FILE="ao12-neighbours.dat", ACTION="READ")
   
   do i=1, npb
    READ(fileID_in3, *) junk1_int, (ao12_neigh(j,i), j=1,12)
   enddo

ELSE
PRINT*, "ao12-neighbours.dat file does not exist"
OPEN(UNIT=fileID_in3, FILE="ao12-neighbours.dat", ACTION="WRITE")

dab = dao_cut
no_neigh=no_ao_neigh
CALL  nearest_neighbour_ab(no_neigh, npb, no, dab, a_atms(:,1:npb), o_atms(:,1:no),ao12_neigh(:,1:npb))
 
!dab = dao_min
!no_neigh=no_ao_neigh
!CALL  nearest_neighbour_ab1(no_neigh, npb, no, dab, a_atms(:,1:npb), o_atms(:,1:no),ao12_neigh(:,1:npb))
 
do i=1, npb
 WRITE(fileID_in3, *) i, (ao12_neigh(j,i), j=1,12)
 enddo

 ENDIF

INQUIRE(FILE="bo6-neighbours.dat", exist=fexist)
IF (fexist) THEN
   PRINT*, "bo6-neighbours.dat file exist"

   OPEN(UNIT=fileID_in4, FILE="bo6-neighbours.dat", ACTION="READ")

   do i=1, nti
    READ(fileID_in4, *) junk1_int, (bo6_neigh(j,i), j=1,6)
   enddo

ELSE
PRINT*, "bo6-neighbours.dat file does not exist"
OPEN(UNIT=fileID_in4, FILE="bo6-neighbours.dat", ACTION="WRITE")


dab = dbo_cut
no_neigh=no_bo_neigh
CALL  nearest_neighbour_ab(no_neigh, nti, no, dab, b_atms(:,1:nti), o_atms(:,1:no), bo6_neigh(:,1:nti)) 

!dab = dbo_min
!no_neigh=no_bo_neigh
!CALL  nearest_neighbour_ab1(no_neigh, nti, no, dab, b_atms(:,1:nti), o_atms(:,1:no), bo6_neigh(:,1:nti)) 


 do i=1, nti
 WRITE(fileID_in4, *) i, (bo6_neigh(j,i), j=1,6)
 enddo

 ENDIF

!print*, 'ok'
CALL CPU_TIME(finish)
PRINT*, "  //////////////// "
PRINT*, " Total CPU time for calculating nearest neighbour and distance =", (finish-start)/60.0, "Minutes"
!
! Calculate a-site polarization
born_charge=zqpb 
CALL cal_a_site_pol_disp(born_charge, npb, no_ao_neigh)



! Calculate b-site polarization  
CALL cal_b_site_pol_disp(zqw, zqmg, nti,  no_bo_neigh)

! For ovito vizualization only 
do i =1, npb 
  do j =1 , 3
  a_atms_symm(j, i) = a_atms(j,i)-a_disp(j,i)
  b_atms_symm(j,i) = b_atms_symm(j,i)-b_disp(j,i)
  enddo
enddo 

!!!!!  AFE sublattice !!!!!
INQUIRE(FILE="afe-sublattice.dat", exist=fexist)
IF (fexist) THEN
   PRINT*, "afe-sublattice.dat file exist"

   OPEN(UNIT=fileID_afe_in, FILE="afe-sublattice.dat", ACTION="READ")

   m1 = nti/2 
   do i=1, m1 
     READ(fileID_afe_in, *) junk1_int, pol_up_index(i), pol_down_index(i)
   enddo
    close(fileID_afe_in)
    nup_spin = m1 
    ndown_spin = m1 
ELSE

  pup(1:3) = a_pol(1:3,1)
  norm = SQRT(pup(1)**2 + pup(2)**2 + pup(3)**2)
  inorm = 1.0/norm
  
  pup(1) = pup(1)*inorm
  pup(2) = pup(2)*inorm
  pup(3) = pup(3)*inorm
  
  m1 =0 
  m2 = 0
  
  DO i =1, npb 
  
  pdown(1:3) = a_pol(1:3,i)
  
  norm = SQRT(pdown(1)**2 + pdown(2)**2 + pdown(3)**2)
  inorm = 1.0/norm
  
  pdown(1) = pdown(1)*inorm
  pdown(2) = pdown(2)*inorm
  pdown(3) = pdown(3)*inorm
  
  angle1 = ACOS(DOT_PRODUCT(pup,pdown))
  
  IF ((angle1-1.5)<0.0) then 
    m1 = m1 + 1
    a_pol_up(1:3, m1) = a_pol(1:3, i)
    pol_up_index(m1) = i 
  else 
    m2 = m2 +1 
    a_pol_up(1:3, m2) = a_pol(1:3, i)
    pol_down_index(m2) = i 
  endif 
  ENDDO 

  nup_spin = m1 
  ndown_spin = m2 
  IF (m1 .NE. m2 ) THEN 
  PRINT*, "UP=", m1, "down = ", m2 
  print*, " Program will stop here"
  STOP 
  else 
    PRINT*, "UP=", m1, "down = ", m2 
ENDIF 
!!!!!  END AFE sublattice !!!!!!!!!!!!!!!!!!!!!!!!!!
!write(fileID_afe, '(A)') ' No   UP-lattice   DOWN-lattice' 

DO i =1 , m1 
  write(fileID_afe, '(1X, 3I5)') i, pol_up_index(i), pol_down_index(i) 
enddo 

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


!CALL cal_afe_order(npb, no_aa_neigh)

! Write file
fileID=12
filename='xcrysden-init.xsf'
CALL write_xsf(fileID, filename,  npb, ntot)
close(12)
!CALL write_pol(npb)
CALL write_pol_afe(nup_spin, ndown_spin, npb)

CALL write_pol_ovito_pmwo(nup_spin, ndown_spin, npb)

!CALL write_dump_ovito( npb)
CALL write_pol_ovito(npb) 

PRINT*, " reading dump >1 "

PRINT*, " $$$$$$$$$ "
PRINT*, " Unit of the Local and Total polarization is in C/m^2"
!PRINT*,  (ptot(j), j = 1, 3), SQRT(ptot(1)**2 + ptot(2)**2 + ptot(3)**2)
CALL CPU_TIME(start)
!!!!!!!! Second itrations!!!!!!!!!!!

DO ndump=2, tot_ndump
! reading dump file
! CALL read_dumpxyz(npb, nti, no)
CALL read_dumpxyz_drew(npb, nti, no)
!  CALL readdump_type1(npb, nti, no)
  ! CALL readdump_type2(npb, nti, no)
   !!*******************************!!
            !  cell volume
    vol = ABS( latt_vec(1,1) * (latt_vec(2,2)*latt_vec(3,3) - latt_vec(2,3)*latt_vec(3,2)) &
               - latt_vec(2,1) * (latt_vec(1,2)*latt_vec(3,3) - latt_vec(1,3)*latt_vec(3,2)) &
               + latt_vec(3,1) * (latt_vec(1,2)*latt_vec(2,3) - latt_vec(1,3)*latt_vec(2,2)) )
    unitvol=vol/nti
   ! Calculate a-site polarization
   born_charge=zqpb 
   CALL cal_a_site_pol_disp(born_charge, npb, no_ao_neigh)
            
  ! Calculate b-site polarization

   CALL cal_b_site_pol_disp(zqw, zqmg, nti,  no_bo_neigh)
    

DO i = 1, nup_spin
  
m1 = pol_up_index(i)
m2 = pol_down_index(i)
a_pol_up(1:3, i) = a_pol(1:3, m1)
a_pol_down(1:3,i) = a_pol(1:3, m2)
ENDDO 
!%%%%%%%%%% 
! pup(1:3) = a_pol(1:3,1)
! norm = SQRT(pup(1)**2 + pup(2)**2 + pup(3)**2)
! inorm = 1.0/norm

! pup(1) = pup(1)*inorm
! pup(2) = pup(2)*inorm
! pup(3) = pup(3)*inorm

! m1 =0 
! m2 = 0
! DO i =1, npb 

! pdown(1:3) = a_pol(1:3,i)

! norm = SQRT(pdown(1)**2 + pdown(2)**2 + pdown(3)**2)
! inorm = 1.0/norm

! pdown(1) = pdown(1)*inorm
! pdown(2) = pdown(2)*inorm
! pdown(3) = pdown(3)*inorm

! angle1 = ACOS(DOT_PRODUCT(pup,pdown))
! IF ((angle1-1.0)<0.0) then 
!   m1 = m1 + 1
!   a_pol_up(1:3, m1) = a_pol(1:3, i)
! else 
!   m2 = m2 +1 
!   a_pol_up(1:3, m2) = a_pol(1:3, i)
! endif 
! ENDDO 

! nup_spin = m1 
! ndown_spin = m2 

!IF (m1 .NE. m2 ) PRINT*, "UP=", m1, "down = ", m2 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  
!CALL cal_afe_order(npb, no_aa_neigh)
 !CALL write_dump_ovito( npb) 

!CALL write_pol(npb)
CALL write_pol_afe(nup_spin, ndown_spin, npb)

CALL write_pol_ovito( npb) 
CALL write_pol_ovito_pmwo(nup_spin, ndown_spin, npb)


    IF ( MOD(ndump,100)==0) THEN
           PRINT*, " number of dumps", ndump, "// in step of 100"
    ENDIF
ENDDO ! ndump

fileID=12
filename='xcrysden-final.xsf'
CALL write_xsf(fileID, filename,  npb, ntot)
close(12)

CALL CPU_TIME(finish)
PRINT*, "  //////////////// "
PRINT*, " Total CPU time =", (finish-start)/3600.0, "hrs"

close(fileID_in)
close(fileID_ovito)
close(fileID_locpol)
close(fileID_totpol)
close(fileID_afe)
close(fileID_disp_a)
close(fileID_disp_b)
close(fileID_rtheta_plt)
END PROGRAM POLARIZATION_CAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
