! PMN_PT
!
PROGRAM polarization_cal
 IMPLICIT NONE

 INTEGER, PARAMETER :: ndim=1000, DP=SELECTED_REAL_KIND(15)
             !!!!! Born effective charge of BiFeO3 #
  !     https://link.aps.org/doi/10.1103/PhysRevB.71.014113
   !!!!! Born effective charge of BaTiO3 #
 
 !https://link.aps.org/doi/10.1103/PhysRevLett.72.3618
! REAL(DP), PARAMETER  :: qpb=2.75d0, qti=7.16d0,&
! qo1=-5.69d0, qo2=-2.11d0, qo3=-2.11d0

 REAL(DP), PARAMETER  :: qa1=4.92d0, qa2=2.75d0, qb1=4.25d0, qb2=7.17d0, &
                          qo1=-3.60d0, qo2=-3.60d0, qo3=-3.60d0

! qa1 = Bi, qa2 = Ba, qb1 = Fe, qb2 = Ti
!                          
 INTEGER :: na, na1, na2, nb1, nb2, no, ntot, ntype, junk, atmnum(5)

 REAL(DP) ::a(3,3), vol, unitvol, p_a1(3), p_a2(3), p_o(3), p_b1(3), p_b2(3), ptot(3), ptemp, &
            atm_a1(ndim,3), atm_a2(ndim,3), atm_b1(ndim,3), atm_b2(ndim,3), o_atm(ndim,3), &
            punit(3), f(3), atm_a(ndim,3) 

 REAL(DP) :: a_disp(ndim,3), a2_disp(ndim,3), b1_disp(ndim,3), b2_disp(ndim,3)
 INTEGER :: b1_a(ndim,8), b2_a(ndim,8)

 REAL :: x, y, z
 INTEGER :: i, j, k, no_neigh
 CHARACTER (LEN=20) :: fname

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fname='struc.vasp'

OPEN(UNIT=12,FILE=fname, ACTION='READ') ! open file for reading
OPEN(UNIT=13,FILE='cell_P.dat',ACTION='WRITE')
OPEN(UNIT=14, FILE='xcryden_P.xsf', ACTION='WRITE')
OPEN(UNIT=15, FILE='total_P.dat', ACTION='WRITE')
OPEN(UNIT=16, FILE='Asite_disp.dat', ACTION='WRITE')
OPEN(UNIT=17, FILE='Bsite_disp.dat', ACTION='WRITE')
!OPEN(UNIT=18, FILE='Mg_disp.dat', ACTION='WRITE')
!OPEN(UNIT=19, FILE='Nb_disp.dat', ACTION='WRITE')
!OPEN(UNIT=20, FILE='Input_struc.xsf', ACTION='WRITE')
!!! Initialization !!!
   DO i=1,3
      DO j=1,3
      a(i,j)=0.0
         ENDDO
             ENDDO

    ntot=0
    na1=0
    na2=0
    nb1=0
    nb2=0
    no=0
    na = 0

    ntype=0
    atmnum(1)=83 !a1,  Bi
    atmnum(2)=56 !a2,  Ba
    atmnum(3)=26 !b1,  Fe
    atmnum(4)=22 !b2,  Ti
    atmnum(5)= 8 ! O
! reading file; file format .xsf
READ(12,*)
READ(12,*)

DO i=1,3
 READ(12,*) a(i,1), a(i,2), a(i,3)
PRINT*, a(i,1), a(i,2), a(i,3)
ENDDO

READ(12,*)
READ(12,*) na1, nb1, nb2, no

READ(12,*)

DO i=1, na1
   READ(12,*) x, y, z
   atm_a1(i,1)=DBLE(x)
   atm_a1(i,2)=DBLE(y)
   atm_a1(i,3)=DBLE(z)
ENDDO
atm_a1 = matmul(atm_a1, a)

!DO i=1, na2
!   READ(12,*) x, y, z
!   atm_a2(i,1)=DBLE(x)
!   atm_a2(i,2)=DBLE(y)
!   atm_a2(i,3)=DBLE(z)
!ENDDO
!atm_a2 = matmul(atm_a2,a)

!!!!!! 
DO i=1, na1
   atm_a(i,1)=atm_a1(i,1)
   atm_a(i,2)=atm_a1(i,2)
   atm_a(i,3)=atm_a1(i,3)
ENDDO
DO i=1, na2
   atm_a(i+na1,1)=atm_a2(i,1)
   atm_a(i+na1,2)=atm_a2(i,2)
   atm_a(i+na1,3)=atm_a2(i,3)
ENDDO
!!!!!
DO i=1, nb1
   READ(12,*) x, y, z
   atm_b1(i,1)=DBLE(x)
   atm_b1(i,2)=DBLE(y)
   atm_b1(i,3)=DBLE(z)
ENDDO
atm_b1 = matmul(atm_b1,a)

DO i=1, nb2
   READ(12,*) x, y, z
   atm_b2(i,1)=DBLE(x)
   atm_b2(i,2)=DBLE(y)
   atm_b2(i,3)=DBLE(z)
ENDDO
atm_b2 = matmul(atm_b2,a)

DO i=1, no
   READ(12,*) x, y, z
   o_atm(i,1)=DBLE(x)
   o_atm(i,2)=DBLE(y)
   o_atm(i,3)=DBLE(z)
ENDDO

o_atm =matmul(o_atm,a)
ntot = na1+nb1+nb2+no  !+na2+nb1+nb2+no
na = na1 !+ na2 
!DO i=1, 7
!READ(12,*)
!ENDDO
!DO i=1,3
! READ(12,*) a(i,1), a(i,2), a(i,3)
! PRINT*, a(i,1), a(i,2), a(i,3)
!
!ENDDO
!READ(12,*)
!READ(12,*) ntot, junk
!
!      IF (ntot>5000) THEN
!          WRITE(*,*) "Please change array length of A, B, C atoms"
!          STOP
!      END IF
!

CLOSE(12)

 !ALLOCATE (pb_disp(npb,3), ti_disp(nti,3), mg_disp(nmg,3), nb_disp(nnb,3),&
!           ti_pb(npb,8), mg_pb(nmg,8), nb_pb(nnb,8))

!  cell volume
          vol = ABS( a(1,1) * (a(2,2)*a(3,3) - a(2,3)*a(3,2)) &
          &         - a(2,1) * (a(1,2)*a(3,3) - a(1,3)*a(3,2)) &
          &         + a(3,1) * (a(1,2)*a(2,3) - a(1,3)*a(2,2)) )
unitvol=vol/(nb1+nb2)
!
PRINT*, na1, nb1, nb2, no

DO i=1, 3
  p_a1(i)=0.0 ! Polarization for A atoms
  p_a2(i)=0.0 ! Polarization for B site
  p_b1(i)=0.0 ! Polarization for B site
  p_b2(i)=0.0 ! Polarization for B site
  p_o(i)=0.0 ! Polarization for O atoms
  f(i)=0.0
ENDDO

!DO i=1, 3
!print*, i
!write(20,*) (a(i,j),j=1,3)
!ENDDO
! Calculating displacement for A toms
no_neigh=12
Call displacement(no_neigh, na1, no, a, atm_a1(1:na1,:), o_atm(1:no,:), a_disp(1:na1,:))
p_a1=(/0.0,0.0,0.0/)

DO i=1, na1
  DO j=1,3
  p_a1(j)=p_a1(j)+ a_disp(i,j)*qa1
   ENDDO
   WRITE(16,200) i, (a_disp(i,j), j=1,3), SQRT(a_disp(i,1)**2 + a_disp(i,2)**2 + a_disp(i,3)**2)
        ENDDO

! Calculating displacement for A toms
!no_neigh=12
!Call displacement(no_neigh, na2, no, a, atm_a2(1:na2,:), o_atm(1:no,:), a2_disp(1:na2,:))
!p_a2=(/0.0,0.0,0.0/)
!
!DO i=1, na2
!  DO j=1,3
!  a_disp(i+na1,j) = a2_disp(i,j)  
!  p_a2(j)=p_a2(j)+ a2_disp(i,j)*qa2
!   ENDDO
!   WRITE(16,200) i, (a2_disp(i,j), j=1,3), SQRT(a2_disp(i,1)**2 + a2_disp(i,2)**2 + a2_disp(i,3)**2)
!        ENDDO        
!
! Calculating displacement for B atoms
no_neigh=6
Call displacement(no_neigh, nb1, no, a, atm_b1(1:nb1,:), o_atm(1:no,:), b1_disp)
!! Polarization due to displacement of ti atom relative to 6 O atoms
p_b1=(/0.0,0.0,0.0/)
DO i=1, nb1
  DO j=1,3
  p_b1(j)=p_b1(j)+ b1_disp(i,j)*qb1
   ENDDO
WRITE(17,200) i, (b1_disp(i,j), j=1,3), SQRT(b1_disp(i,1)**2 + b1_disp(i,2)**2 + b1_disp(i,3)**2)
        ENDDO

no_neigh=6
Call displacement(no_neigh, nb2, no, a, atm_b2(1:nb2,:), o_atm(1:no,:), b2_disp)
!! Polarization due to displacement of ti atom relative to 6 O atoms
p_b2=(/0.0,0.0,0.0/)
DO i=1, nb2
  DO j=1,3
  p_b2(j)=p_b2(j)+ b2_disp(i,j)*qb2
   ENDDO
WRITE(17,200) i, (b2_disp(i,j), j=1,3), SQRT(b2_disp(i,1)**2 + b2_disp(i,2)**2 + b2_disp(i,3)**2)
ENDDO


200 FORMAT(I3, 4F12.5)
!
!!      (1.602*10^(-19)c/Angs^2)=16.02* c/m^2 unit conversion
DO i=1, 3
!PRINT*, p_pb(i), p_ti(i), p_mg(i), p_nb(i)

 !ptot(i)= ((p_a1(i)+p_a2(i) + p_b1(i) + p_b2(i))*16.02)/vol  ! c/m^2
 ptot(i) = ((p_a1(i)+p_b1(i)+p_b2(i))*16.02)/vol  ! c/m^2
 ENDDO
!
WRITE(15,*) '#total polarization along x, y, z in unit of C/m^2'
WRITE(15, 400) ptot(1), ptot(2), ptot(3), SQRT(ptot(1)**2 + ptot(2)**2 + ptot(3)**2)
400 FORMAT (1X, 4F12.8)
PRINT*, 'Total Polarization along x, y, z in unit of C/m^2'
PRINT*, ptot(1), ptot(2), ptot(3), SQRT(ptot(1)**2 + ptot(2)**2 + ptot(3)**2)
!
!calculating number of nearest-neighbour A atoms around a B in a unit cell
no_neigh=8
CALL  nearest_neighbour(no_neigh, nb1, na, a, atm_b1(1:nb1,:), atm_a(1:na,:),b1_a)

no_neigh=8
CALL  nearest_neighbour(no_neigh, nb2, na, a, atm_b2(1:nb2,:), atm_a(1:na,:),b2_a)

!!!!!!!WRITEING for xcruden .xsf file FORMAT
WRITE(14,*) '#Pb Ti O'
WRITE(14,*)'CRYSTAL'
WRITE(14,*)'PRIMVEC'
DO i=1,3
  WRITE(14,300)a(i,1), a(i,2), a(i,3)
ENDDO
WRITE(14,*)'CONVVEC'
DO i=1,3
  WRITE(14,300)a(i,1), a(i,2), a(i,3)
ENDDO
300 FORMAT(1X, 3F13.8)
WRITE(14,*)'PRIMCOORD'
WRITE(14,*) ntot, '1'
DO i=1, na1
f(1)=a_disp(i,1)*qa1/unitvol*16.02
f(2)=a_disp(i,2)*qa1/unitvol*16.02
f(3)=a_disp(i,3)*qa1/unitvol*16.02

  WRITE(14,301)atmnum(1), atm_a(i,1), atm_a(i,2), atm_a(i,3), f(1), f(2),f(3)
ENDDO

!DO i=1, na2
!f(1)=a_disp(i+na1,1)*qa2/unitvol*16.02
!f(2)=a_disp(i+na1,2)*qa2/unitvol*16.02
!f(3)=a_disp(i+na1,3)*qa2/unitvol*16.02
!
!  WRITE(14,301)atmnum(2), atm_a(i+na1,1), atm_a(i+na1,2), atm_a(i+na1,3), f(1), f(2),f(3)
!ENDDO
!
!!Polarization per unit cell
DO i=1, nb1
  DO j=1,3 ! x y z
  ptemp=0.0
    DO k=1,8
      ptemp=ptemp+ a_disp(b1_a(i,k),j)*qa1  ! B atom Polarization
    ENDDO
!    PRINT*, ptemp/8.0, pb(j), pa(j)
 punit(j)=(b1_disp(i,j)*qb1 + (ptemp/8.0))/unitvol*16.02  ! Total polarization per unit cell
ENDDO ! x y z

WRITE(13,302) (atm_b1(i,j),j=1,3), (punit(j),j=1,3)

f(1)=b1_disp(i,1)*qb1/unitvol*16.02
f(2)=b1_disp(i,2)*qb1/unitvol*16.02
f(3)=b1_disp(i,3)*qb1/unitvol*16.02

WRITE(14,301) atmnum(3), (atm_b1(i,j),j=1,3), f(1), f(2), f(3) !(punit(j),j=1,3)
ENDDO

DO i=1, nb2
  DO j=1,3 ! x y z
  ptemp=0.0
    DO k=1,8
      ptemp=ptemp+ a_disp(b2_a(i,k),j)*qa2  ! B atom Polarization
    ENDDO
!    PRINT*, ptemp/8.0, pb(j), pa(j)
 punit(j)=(b2_disp(i,j)*qb2 + (ptemp/8.0))/unitvol*16.02  ! Total polarization per unit cell
ENDDO ! x y z

WRITE(13,302) (atm_b2(i,j),j=1,3), (punit(j),j=1,3)

f(1)=b2_disp(i,1)*qb2/unitvol*16.02
f(2)=b2_disp(i,2)*qb2/unitvol*16.02
f(3)=b2_disp(i,3)*qb2/unitvol*16.02

WRITE(14,301) atmnum(4), (atm_b2(i,j),j=1,3), f(1), f(2), f(3) !(punit(j),j=1,3)
ENDDO

f(1)=0.0
f(2)=0.0
f(3)=0.0
!
DO i=1, no
  WRITE(14,301) atmnum(5), o_atm(i,1), o_atm(i,2), o_atm(i,3), f(1), f(2), f(3)
ENDDO
301 FORMAT (1X, I5, 6F13.8)
302 FORMAT (1X,  6F13.8)
       END PROGRAM polarization_cal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE displacement(n1st, n1, n2, latt, atm1, atm2, disp)
        IMPLICIT NONE
        INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(15)
        INTEGER, INTENT(IN) :: n1st, n1, n2
        REAL(dp), INTENT(IN) :: latt(3,3), atm1(n1,3), atm2(n2,3)
        REAL(dp), INTENT(OUT) :: disp(n1, 3)

        REAL(dp) :: v1(3), vi(3), vf(3), atm_info(27*n2,4), d(27*n2),dmax, &
                     c_mass(3), neigh(n1st, 3)
!
        INTEGER :: i, j, k, n, l1, l2,l3, loc

         dmax=100 ! angs

        DO i=1, n1
            v1(1)=atm1(i,1)
            v1(2)=atm1(i,2)
            v1(3)=atm1(i,3)
            n=0
            DO j=1, n2
               vi(1)=atm2(j,1)
               vi(2)=atm2(j,2)
               vi(3)=atm2(j,3)
               DO l1=-1, 1
                  DO l2 =-1, 1
                     DO l3=-1, 1
                       n=n+1
                        vf(1)=vi(1)+ l1*latt(1,1) + l2*latt(2,1) + l3*latt(3,1)!
                        vf(2)=vi(2)+ l2*latt(2,2) + l1*latt(1,2) + l3*latt(3,2)!
                        vf(3)=vi(3)+ l3*latt(3,3) + l1*latt(1,3) + l2*latt(2,3)

                        atm_info(n,1)=j
                        atm_info(n,2)=vf(1)
                        atm_info(n,3)=vf(2)
                        atm_info(n,4)=vf(3)
                        d(n)=SQRT((v1(1)-vf(1))**2 + &
                                  (v1(2)-vf(2))**2 + &
                                  (v1(3)-vf(3))**2)
                       ENDDO
                     ENDDO
                   ENDDO
               ENDDO
!
                 DO k=1, n1st
                    loc=MINLOC(d, DIM=1)
                     neigh(k,1)=atm_info(loc,2)
                     neigh(k,2)=atm_info(loc,3)
                     neigh(k,3)=atm_info(loc,4)
                     !ab_loc(i,k)=INT(atm_info(loc,1))
                     d(loc)=dmax ! to find new min value location
                     ENDDO
        ! Center of mass
                   c_mass(1)=SUM(neigh(1:n1st,1))/DBLE(n1st)
                   c_mass(2)=SUM(neigh(1:n1st,2))/DBLE(n1st)
                   c_mass(3)=SUM(neigh(1:n1st,3))/DBLE(n1st)
       ! displacement
              disp(i,1)=v1(1)-c_mass(1)
              disp(i,2)=v1(2)-c_mass(2)
              disp(i,3)=v1(3)-c_mass(3)

                    ENDDO
END SUBROUTINE displacement

SUBROUTINE nearest_neighbour(n1st, n1, n2, latt, atm1, atm2, ab_loc)
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(15)
  INTEGER, INTENT(IN) :: n1st, n1, n2
  REAL(dp), INTENT(IN) :: latt(3,3), atm1(n1,3), atm2(n2,3)
  INTEGER, INTENT(OUT) :: ab_loc(n1, n1st)

  REAL(dp) :: v1(3), vi(3), vf(3), atm_info(27*n2,4), d(27*n2),dmax
!
  INTEGER :: i, j, k, n, l1, l2, l3, loc

   dmax=100 ! angs

  DO i=1, n1
      v1(1)=atm1(i,1)
      v1(2)=atm1(i,2)
      v1(3)=atm1(i,3)
      n=0
      DO j=1, n2
         vi(1)=atm2(j,1)
         vi(2)=atm2(j,2)
         vi(3)=atm2(j,3)
         DO l1=-1, 1
            DO l2 =-1, 1
               DO l3=-1, 1
                 n=n+1
                 vf(1)=vi(1)+ l1*latt(1,1) + l2*latt(2,1) + l3*latt(3,1)!
                 vf(2)=vi(2)+ l2*latt(2,2) + l1*latt(1,2) + l3*latt(3,2)!
                 vf(3)=vi(3)+ l3*latt(3,3) + l1*latt(1,3) + l2*latt(2,3)

                  atm_info(n,1)=j
                  atm_info(n,2)=vf(1)
                  atm_info(n,3)=vf(2)
                  atm_info(n,4)=vf(3)
                  d(n)=SQRT((v1(1)-vf(1))**2 + &
                            (v1(2)-vf(2))**2 + &
                            (v1(3)-vf(3))**2)
                 ENDDO
               ENDDO
             ENDDO
         ENDDO
!
           DO k=1, n1st
              loc=MINLOC(d, DIM=1)
               ab_loc(i,k)=INT(atm_info(loc,1))
               d(loc)=dmax ! to find new min value location
               ENDDO
          ENDDO

END SUBROUTINE nearest_neighbour
