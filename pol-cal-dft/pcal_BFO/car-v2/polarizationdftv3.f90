! PMN_PT
!
     PROGRAM polarization_cal
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(15), dp_int=SELECTED_INT_KIND(7)
      INTEGER(dp_int), PARAMETER :: ndim = 100000
                  !!!!! Born effective charge of BaTiO3 #
             !https://link.aps.org/doi/10.1103/PhysRevLett.72.3618
      REAL(DP), PARAMETER  :: qpb=2.75d0, qti=7.16d0, qmg=2.0, qnb=7.5, &
                               qo1=-5.69d0, qo2=-2.11d0, qo3=-2.11d0


      INTEGER :: npb, nti, nmg, nnb, no, ntot, ntype, junk, atmnum(5)

      REAL(DP) ::a(3,3), vol, unitvol, p_pb(3), p_ti(3), p_o(3), p_mg(3), p_nb(3), ptot(3), ptemp, &
                 pb_atm(ndim,3), ti_atm(ndim,3), mg_atm(ndim,3), nb_atm(ndim,3), o_atm(ndim,3), &
               punit(3),f(3)

      REAL(DP), ALLOCATABLE :: pb_disp(:,:), ti_disp(:,:), mg_disp(:,:), nb_disp(:,:)
      INTEGER, ALLOCATABLE :: ti_pb(:,:), mg_pb(:,:), nb_pb(:,:)

      REAL :: x, y, z
      INTEGER :: i, j, k, no_neigh
      CHARACTER (LEN=20) :: fname

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fname='struc.vasp'

OPEN(UNIT=12,FILE=fname, ACTION='READ') ! open file for reading
OPEN(UNIT=13,FILE='cell_P.dat',ACTION='WRITE')
OPEN(UNIT=14, FILE='xcryden_P.xsf', ACTION='WRITE')
OPEN(UNIT=15, FILE='total_P.dat', ACTION='WRITE')
OPEN(UNIT=16, FILE='Pb_disp.dat', ACTION='WRITE')
OPEN(UNIT=17, FILE='Ti_disp.dat', ACTION='WRITE')
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
    npb=0
    nti=0
    nmg=0
    nnb=0
    no=0
    ntype=0
    atmnum(1)=82 ! Pb
    atmnum(2)=22 ! Ti
    atmnum(3)=8  ! O
    atmnum(4)=12 !Mg
    atmnum(5)=41 !Nb
! reading file; file format .xsf
READ(12,*)
READ(12,*)

DO i=1,3
 READ(12,*) a(i,1), a(i,2), a(i,3)
PRINT*, a(i,1), a(i,2), a(i,3)
ENDDO

READ(12,*)
READ(12,*) npb, nti, no

READ(12,*)

DO i=1, npb
   READ(12,*) x, y, z
   pb_atm(i,1)=DBLE(x)
   pb_atm(i,2)=DBLE(y)
   pb_atm(i,3)=DBLE(z)
ENDDO

DO i=1, nti
   READ(12,*) x, y, z
   ti_atm(i,1)=DBLE(x)
   ti_atm(i,2)=DBLE(y)
   ti_atm(i,3)=DBLE(z)
ENDDO

DO i=1, no
   READ(12,*) x, y, z
   o_atm(i,1)=DBLE(x)
   o_atm(i,2)=DBLE(y)
   o_atm(i,3)=DBLE(z)
ENDDO

ntot=npb+nti+no

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
!
!DO i=1, ntot
!  READ(12,*) ntype, x, y, z
! IF (ntype==atmnum(1)) THEN ! Pb atoms
!    npb=npb+1
!      pb_atm(npb,1)=DBLE(x)
!      pb_atm(npb,2)=DBLE(y)
!      pb_atm(npb,3)=DBLE(z)
!  ELSEIF(ntype==atmnum(2)) THEN ! Ti atoms
!     nti=nti+1
!     ti_atm(nti,1)=DBLE(x)
!     ti_atm(nti,2)=DBLE(y)
!     ti_atm(nti,3)=DBLE(z)
!   ELSEIF(ntype==atmnum(3))THEN ! O atoms
!     no=no+1
!      o_atm(no,1)=DBLE(x)
!      o_atm(no,2)=DBLE(y)
!      o_atm(no,3)=DBLE(z)
!   ELSEIF(ntype==atmnum(4)) THEN ! Mg atoms
!     nmg=nmg+1
!     mg_atm(nmg,1)=DBLE(x)
!     mg_atm(nmg,2)=DBLE(y)
!     mg_atm(nmg,3)=DBLE(z)
!  ELSEIF(ntype==atmnum(5)) THEN ! Nb atoms
!     nnb=nnb+1
!     nb_atm(nnb,1)=DBLE(x)
!     nb_atm(nnb,2)=DBLE(y)
!     nb_atm(nnb,3)=DBLE(z)
!ELSE
!      PRINT*, 'Error:: No species more than three type'
!      PRINT*, " Please change code"
!      STOP
!  ENDIF
!
! ENDDO
CLOSE(12)

 ALLOCATE (pb_disp(npb,3), ti_disp(nti,3), mg_disp(nmg,3), nb_disp(nnb,3),&
           ti_pb(npb,8), mg_pb(nmg,8), nb_pb(nnb,8))

!  cell volume
          vol = ABS( a(1,1) * (a(2,2)*a(3,3) - a(2,3)*a(3,2)) &
          &         - a(2,1) * (a(1,2)*a(3,3) - a(1,3)*a(3,2)) &
          &         + a(3,1) * (a(1,2)*a(2,3) - a(1,3)*a(2,2)) )
unitvol=vol/(nti+nmg+nnb)
!
PRINT*, npb, nti, no

DO i=1, 3
  p_pb(i)=0.0 ! Polarization for A atoms
  p_ti(i)=0.0 ! Polarization for B site
  p_mg(i)=0.0 ! Polarization for B site
  p_nb(i)=0.0 ! Polarization for B site
  p_o(i)=0.0 ! Polarization for O atoms
  f(i)=0.0
ENDDO

!DO i=1, 3
!print*, i
!write(20,*) (a(i,j),j=1,3)
!ENDDO
! Calculating displacement for A toms
no_neigh=12
Call displacement(no_neigh, npb, no, a, pb_atm(1:npb,:), o_atm(1:no,:), pb_disp)
p_pb=(/0.0,0.0,0.0/)

DO i=1, npb
  DO j=1,3
  p_pb(j)=p_pb(j)+ pb_disp(i,j)*qpb
   ENDDO
   WRITE(16,200) i, (pb_disp(i,j), j=1,3), SQRT(pb_disp(i,1)**2 + pb_disp(i,2)**2 + pb_disp(i,3)**2)
        ENDDO

! Calculating displacement for B atoms
no_neigh=6
Call displacement(no_neigh, nti, no, a, ti_atm(1:nti,:), o_atm(1:no,:), ti_disp)
!! Polarization due to displacement of ti atom relative to 6 O atoms
p_ti=(/0.0,0.0,0.0/)
DO i=1, nti
  DO j=1,3
  p_ti(j)=p_ti(j)+ ti_disp(i,j)*qti
   ENDDO
WRITE(17,200) i, (ti_disp(i,j), j=1,3), SQRT(ti_disp(i,1)**2 + ti_disp(i,2)**2 + ti_disp(i,3)**2)
        ENDDO
!
!!! Calculating displacement for B prime atoms
!no_neigh=6
!Call displacement(no_neigh, nmg, no, a, mg_atm(1:nmg,:), o_atm(1:no,:), mg_disp)
!!!! Polarization due to displacement of Pb atom relative to 6 O atoms
!p_mg=(/0.0,0.0,0.0/)
!DO i=1, nmg
!  DO j=1,3
!  p_mg(j)=p_mg(j)+ mg_disp(i,j)*qmg
!   ENDDO
!WRITE(18,200) i, (mg_disp(i,j), j=1,3), SQRT(mg_disp(i,1)**2 + mg_disp(i,2)**2 + mg_disp(i,3)**2)
!        ENDDO
!!
!!! Calculating displacement for B prime atoms
!no_neigh=6
!Call displacement(no_neigh, nnb, no, a, nb_atm(1:nnb,:), o_atm(1:no,:), nb_disp)
!!! Polarization due to displacement of Pb atom relative to 6 O atoms
!p_nb=(/0.0,0.0,0.0/)
!DO i=1, nnb
!  DO j=1,3
!  p_nb(j)=p_nb(j)+ nb_disp(i,j)*qnb
!
!   ENDDO
!WRITE(19,200) i, (nb_disp(i,j), j=1,3), SQRT(nb_disp(i,1)**2 + nb_disp(i,2)**2 + nb_disp(i,3)**2)
!        ENDDO
!!
200 FORMAT(I3, 4F12.5)
!
!!      (1.602*10^(-19)c/Angs^2)=16.02* c/m^2 unit conversion
DO i=1, 3
!PRINT*, p_pb(i), p_ti(i), p_mg(i), p_nb(i)

 ptot(i)= ((p_pb(i)+p_ti(i))*16.02)/vol  ! c/m^2 
 !((p_pb(i)+p_ti(i)+p_mg(i)+p_nb(i) )*16.02)/vol  ! c/m^2
 ENDDO
!
WRITE(15,*) '#total polarization along x, y, z in unit of C/m^2'
WRITE(15, 400) ptot(1), ptot(2), ptot(3)
400 FORMAT (1X, 3F12.8)
PRINT*, 'Total Polarization along x, y, z in unit of C/m^2'
PRINT*, 'p1=', ptot(1),'p2=', ptot(2), 'p3=', ptot(3), '|P|=', SQRT(ptot(1)**2 + ptot(2)**2 + ptot(3)**2)
!
!calculating number of nearest-neighbour A atoms around a B in a unit cell
no_neigh=8
CALL  nearest_neighbour(no_neigh, nti, npb, a, ti_atm(1:nti,:), pb_atm(1:npb,:), ti_pb)
!CALL  nearest_neighbour(no_neigh, nmg, npb, a, mg_atm(1:nmg,:), pb_atm(1:npb,:), mg_pb)
!CALL  nearest_neighbour(no_neigh, nnb, npb, a, nb_atm(1:nnb,:), pb_atm(1:npb,:), nb_pb)
!!DO i=1, nb
!!  PRINT*, (ab(i,j),j=1,8)
!!ENDDO
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
DO i=1, npb
  WRITE(14,301)atmnum(1), pb_atm(i,1), pb_atm(i,2), pb_atm(i,3), f(1), f(2), f(3)
ENDDO
!!Polarization per unit cell
DO i=1, nti
  DO j=1,3 ! x y z
  ptemp=0.0
    DO k=1,8
      ptemp=ptemp+ pb_disp(ti_pb(i,k),j)*qpb  ! B atom Polarization
    ENDDO
!    PRINT*, ptemp/8.0, pb(j), pa(j)
 punit(j)=(ti_disp(i,j)*qti + (ptemp/8.0))/unitvol*16.02  ! Total polarization per unit cell
ENDDO ! x y z

WRITE(13,302) (ti_atm(i,j),j=1,3), (punit(j),j=1,3)
WRITE(14,301) atmnum(2), (ti_atm(i,j),j=1,3), (punit(j),j=1,3)
ENDDO

!DO i=1, nmg
!  DO j=1,3
!  ptemp=0.0
!    DO k=1,8
!      ptemp=ptemp+ pb_disp(mg_pb(i,k),j)*qpb  ! B atom Polarization
!    ENDDO
!!    PRINT*, ptemp/8.0, pb(j), pa(j)
! punit(j)=(mg_disp(i,j)*qmg + (ptemp/8.0))/unitvol*16.02  ! Total polarization per unit cell
!ENDDO
!WRITE(13,302) (mg_atm(i,j),j=1,3), (punit(j),j=1,3)
!WRITE(14,301) atmnum(4), (mg_atm(i,j),j=1,3), (punit(j),j=1,3)
!ENDDO
!
!DO i=1, nnb
!  DO j=1,3
!  ptemp=0.0
!    DO k=1,8
!      ptemp=ptemp+ pb_disp(nb_pb(i,k),j)*qpb  ! B atom Polarization
!    ENDDO
!!!    PRINT*, ptemp/8.0, pb(j), pa(j)
! punit(j)=(nb_disp(i,j)*qnb + (ptemp/8.0))/unitvol*16.02  ! Total polarization per unit cell
!ENDDO
!WRITE(13,302) (nb_atm(i,j),j=1,3), (punit(j),j=1,3)
!WRITE(14,301) atmnum(5), (nb_atm(i,j),j=1,3), (punit(j),j=1,3)
!ENDDO
!
DO i=1, no
  WRITE(14,301) atmnum(3), o_atm(i,1), o_atm(i,2), o_atm(i,3), f(1), f(2), f(3)
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
                        vf(1)=vi(1)+ l1*latt(1,1) !
                        vf(2)=vi(2)+ l2*latt(2,2) !
                        vf(3)=vi(3)+ l3*latt(3,3)

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
                  vf(1)=vi(1)+ l1*latt(1,1) !
                  vf(2)=vi(2)+ l2*latt(2,2) !
                  vf(3)=vi(3)+ l3*latt(3,3)

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
