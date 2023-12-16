
 MODULE my_constants
 IMPLICIT NONE
 INTEGER, PARAMETER :: dp_real=SELECTED_REAL_KIND(16), dp_int=SELECTED_INT_KIND(15)

 REAL(dp_real), PARAMETER :: charge_A=1.93753, charge_B=1.77706, charge_O=-1.2382 ! A=Bi, B=Fe, O=O BiFeO3 compound
! REAL(dp_real), PARAMETER :: charge_A=1.34730, charge_B=1.28905, charge_O=-0.87878 ! A=Ba, B=Ti, O=O BaTiO3 compound
! REAL(dp_real), PARAMETER :: charge_A=1.38177, charge_B=0.99997, charge_O=-0.79391 ! A=Pb, B=Ti, O=O PbTiO3 compound
 INTEGER(dp_int), PARAMETER :: nx=1*5,ny=1*5,nz=1*5 ! cell dimension

 REAL, PARAMETER :: massA=208.9804, massB=55.8450,massO=15.9994 !massA=Bi, massB=Fe, massO=O
 !REAL, PARAMETER :: massA=137.327, massB=47.88,massO=15.9994 !massA=Ba, massB=Ti, massO=O
 !REAL, PARAMETER :: massA=207.2, massB=47.88,massO=15.9994 ! massA=Pb, massB=Ti, massO=O
 !INTEGER, PARAMETER :: angle_coeffs=200  !BaTiO3
 !INTEGER, PARAMETER :: angle_coeffs=50  !PaTiO3
 REAL(dp_real), PARAMETER :: strainx=0.0,strainy=0.0, strainz=0.0    !strain in percentage
 !REAL(dp_real):: a11=4.065075, a22=4.065075, a33=4.109  ! a=b=c=4.005 BaTiO3 latice parameters
 REAL(dp_real):: a11=7.782935000, a22=7.782935000, a33=7.782935000  ! a=b=c=3.86 PbTiO3 latice parameters

INTEGER(dp_int), PARAMETER :: natms_max=100000000
INTEGER, PARAMETER :: ndim3=3
INTEGER, PARAMETER :: ndim4=4

END MODULE my_constants

MODULE data_structure
USE my_constants
REAL(dp_real),dimension(:, :), ALLOCATABLE :: a_atms,b_atms,&
o1_atms, o2_atms, o3_atms

REAL(dp_real), dimension(:, :), ALLOCATABLE :: a_atms_opt,b_atms_opt,  &
o1_atms_opt, o2_atms_opt, o3_atms_opt

INTEGER(dp_int), dimension(:, :),ALLOCATABLE :: angle_ti_o_ti
INTEGER(dp_int), dimension(:, :),ALLOCATABLE :: all_angle_ti_o_ti
INTEGER(dp_int), dimension(:),ALLOCATABLE :: fe_bfo

REAL(dp_real), dimension(ndim3, ndim3) :: a_unit, latt_sup
END MODULE data_structure

MODULE data_cubic40
USE my_constants
USE data_structure
IMPLICIT NONE
private
public :: supercell_cubic, neighbour_atm1_atm2,neighbour_axis,&
         projection_car_axis,dm_Zint_angle

CONTAINS

SUBROUTINE supercell_cubic(nu, nuti, atmu, flag, ntotal)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nu, nuti, flag
    REAL (dp_real), INTENT(IN) :: atmu(ndim3, nu)
    INTEGER(dp_int), INTENT(OUT) :: ntotal

    REAL(dp_real), dimension(ndim3, nuti) :: a1_atmu, b1_atmu, &
                      o11_atmu, o22_atmu, o33_atmu

    REAL(dp_real) :: T(3)

    REAL(dp_real) :: a, b, c

    INTEGER(dp_int)  :: ix, iy, iz, n
    INTEGER :: i 

    n = 0
    DO i =1, nuti
       n = n + 1
    a1_atmu(1:3,i) = atmu(1:3, n)
    ENDDO

    DO i =1, nuti
        n = n + 1
        b1_atmu(1:3,i) = atmu(1:3, n)
    ENDDO

    DO i = 1, nuti
     n = n + 1
    o11_atmu(1:3,i) = atmu(1:3, n)
    ENDDO

    DO i = 1, nuti
      n = n + 1
       o22_atmu(1:3,i) = atmu(1:3, n)
    ENDDO

    DO i=1, nuti
    n = n + 1
    o33_atmu(1:3,i) = atmu(1:3, n)
    ENDDO

    a  = a_unit(1,1)
    b  = a_unit(2,2)
    c  = a_unit(3,3)
    n = 0

    ntotal = nx*ny*nz*nuti 

   DO i = 1, nuti

    DO iz =0, nz-1
        DO iy = 0, ny-1
            DO ix = 0, nx-1

            T=(/DBLE(ix)*a, DBLE(iy)*b, DBLE(iz)*c/)     
                n = n + 1    
            if (flag==1) then 
            a_atms(1,n) = DBLE(ix)
            b_atms(1,n) = DBLE(ix)
            o1_atms(1,n) = DBLE(ix)
            o2_atms(1,n) = DBLE(ix)
            o3_atms(1,n) = DBLE(ix)

            a_atms(2,n) = a1_atmu(1,i)*a + T(1)
            a_atms(3,n) = a1_atmu(2,i)*b + T(2)
            a_atms(4,n) = a1_atmu(3,i)*c + T(3)
            

            b_atms(2,n) = b1_atmu(1,i)*a + T(1)
            b_atms(3,n) = b1_atmu(2,i)*b + T(2)
            b_atms(4,n) = b1_atmu(3,i)*c + T(3)

            !print*, b_atms(n,2), b_atms(n, 3), b_atms(n, 4), a, b, c,  b_atmu(i,1),  b_atmu(i,2), b_atmu(i,3)

            o1_atms(2,n) = o11_atmu(1,i)*a + T(1)
            o1_atms(3,n) = o11_atmu(2,i)*b + T(2)
            o1_atms(4,n) = o11_atmu(3,i)*c + T(3)


            o2_atms(2,n) = o22_atmu(1,i)*a + T(1)
            o2_atms(3,n) = o22_atmu(2,i)*b + T(2)
            o2_atms(4,n) = o22_atmu(3,i)*c + T(3)


            o3_atms(2,n) = o33_atmu(1,i)*a  + T(1)
            o3_atms(3,n) = o33_atmu(2,i)*b  + T(2)
            o3_atms(4,n) = o33_atmu(3,i)*c  + T(3)

            elseif (flag==2)then

                 a_atms_opt(1,n) = DBLE(ix)
                 b_atms_opt(1,n) = DBLE(ix)
                o1_atms_opt(1,n) = DBLE(ix)
                o2_atms_opt(1,n) = DBLE(ix)
                o3_atms_opt(1,n) = DBLE(ix)
    
                a_atms_opt(2,n) = a1_atmu(1,i)*a + T(1)
                a_atms_opt(3,n) = a1_atmu(2,i)*b + T(2)
                a_atms_opt(4,n) = a1_atmu(3,i)*c + T(3)
                
    
                b_atms_opt(2,n) = b1_atmu(1,i)*a + T(1)
                b_atms_opt(3,n) = b1_atmu(2,i)*b + T(2)
                b_atms_opt(4,n) = b1_atmu(3,i)*c + T(3)
    
                !print*, b_atms(n,2), b_atms(n, 3), b_atms(n, 4), a, b, c,  b_atmu(i,1),  b_atmu(i,2), b_atmu(i,3)
    
                o1_atms_opt(2,n) = o11_atmu(1,i)*a + T(1)
                o1_atms_opt(3,n) = o11_atmu(2,i)*b + T(2)
                o1_atms_opt(4,n) = o11_atmu(3,i)*c + T(3)
    
    
                o2_atms_opt(2,n) = o22_atmu(1,i)*a + T(1)
                o2_atms_opt(3,n) = o22_atmu(2,i)*b + T(2)
                o2_atms_opt(4,n) = o22_atmu(3,i)*c + T(3)
    
    
                o3_atms_opt(2,n) = o33_atmu(1,i)*a  + T(1)
                o3_atms_opt(3,n) = o33_atmu(2,i)*b  + T(2)
                o3_atms_opt(4,n) = o33_atmu(3,i)*c  + T(3)
            else 
                print*, "===Error==="
            endif

                ENDDO
            ENDDO
        ENDDO
     ENDDO
END SUBROUTINE supercell_cubic



SUBROUTINE neighbour_atm1_atm2(dmin, nsti, atm1, atm2, neigh_ab)
IMPLICIT NONE 
INTEGER(dp_int), INTENT(IN) :: nsti
REAL (dp_real) :: dmin
REAL(dp_real), dimension(ndim4,nsti), INTENT(IN) ::  atm1, atm2
INTEGER(dp_int), INTENT(OUT) :: neigh_ab(8,nsti)
REAL(dp_real) :: v1(3), v2(3), dist
INTEGER(dp_int) :: i, j

INTEGER :: m


DO i =1, nsti
    DO j = 1, 8
        neigh_ab(j,i) = -1
    ENDDO
ENDDO

DO i =1, nsti 
    v1= atm1(2:4,i)
    m =0

    DO j =1, nsti 
        v2 = atm2(2:4,j)

        dist = SQRT(( v1(1) - v2(1))**2 + &
               (v1(2) - v2(2))**2   + &
               (v1(3) - v2(3))**2)

        IF (dist < 1.0) THEN
            CYCLE
        ENDIF  

        if ( (dist-dmin) < 0.0) THEN 
            !print*, dist
            m = m +1
            neigh_ab(m,i) = j
            ENDIF
    ENDDO 
ENDDO

END SUBROUTINE neighbour_atm1_atm2


SUBROUTINE neighbour_axis(vec_unit, dmin, latt_sup, nsti, atm1, atm2, nij, atmsij, neighij)
    IMPLICIT NONE 
    INTEGER(dp_int), INTENT(IN) :: nsti, nij
    REAL(dp_real), INTENT(IN) :: dmin, atm1(4, nsti), atm2(4, nsti)
    REAL(dp_real), INTENT(IN) :: latt_sup(3,3)
    REAL (dp_real), INTENT(IN) :: vec_unit(3)
    INTEGER(dp_int), INTENT(IN) :: atmsij(nsti)
    INTEGER(dp_int), INTENT(OUT) :: neighij(nsti)

REAL(dp_real) :: v1(3), v2(3), w(3), dist, T(3), sl

INTEGER(dp_int) :: i, j, m
INTEGER ::  k



DO i =1, nsti
    neighij(i) = 0
ENDDO

DO i =1, nij
    m = atmsij(i)
    v1= atm1(2:4,m)

    DO j =1, nsti 
        v2 = atm2(2:4,j)

        DO k = -1, 1
           IF (k==0) CYCLE 
            !print*,k
           sl = DBLE(k)
           T(1) =  vec_unit(1)*latt_sup(1,1)*sl
           T(2) =  vec_unit(2)*latt_sup(2,2)*sl
           T(3)  = vec_unit(3)*latt_sup(3,3)*sl

           w(1) = v2(1) + T(1)
           w(2) = v2(2) + T(2)
           w(3) = v2(3) + T(3)

           dist = SQRT((v1(1)-w(1))**2 + (v1(2)-w(2))**2 + (v1(3)-w(3))**2)
           !print*, dist 
           if ( (dist-dmin) < 0.0) THEN 
            !print*, dist
            neighij(m) = j
                ENDIF
            ENDDO

    ENDDO 
ENDDO
END SUBROUTINE neighbour_axis

SUBROUTINE projection_car_axis(vec_unit, nsti, atm1, atm2, neigh_ab, NeighAlongAxis)
IMPLICIT NONE 

INTEGER(dp_int), INTENT(IN) :: nsti, neigh_ab(8,nsti)

REAL(dp_real), dimension(ndim4, nsti), INTENT(IN) :: atm1, atm2
REAL(dp_real), INTENT(IN) :: vec_unit(3)
INTEGER(dp_int), INTENT(OUT) :: NeighAlongAxis(nsti)

REAL(dp_real) :: v1(3), v2(3),  w(3), d, d_proj

INTEGER(dp_int) :: i, j, m, ncout  

DO i = 1, nsti
NeighAlongAxis(i) = 0
ENDDO

DO i = 1, 3
v1(i) = 0.0
v2(i) = 0.0
w(i) = 0.0
ENDDO

DO i =1, nsti 
v1 = atm1(2:4, i)
ncout = 0
DO j = 1, 8
  m = neigh_ab(j, i)
  IF(m==-1) CYCLE
    v2 = atm2(2:4,m)

    w(1) = v2(1)-v1(1)
    w(2) = v2(2)-v1(2)
    w(3) = v2(3)-v1(3)

    d = SQRT(w(1)**2 + w(2)**2 + w(3)**2)

    w(1) = w(1)/d
    w(2) = w(2)/d
    w(3) =  w(3)/d

   d_proj = DOT_PRODUCT(w,vec_unit)

  IF (d_proj > 0.5) THEN
    ncout = ncout + 1
    NeighAlongAxis(i) = m
  ENDIF
ENDDO

IF (ncout>1) PRINT*, "Wrong counting"
ENDDO


END SUBROUTINE projection_car_axis


SUBROUTINE dm_Zint_angle(nsti, vec_unit, latt_vec, b11_atms, o11_atms, angle1_ti_o_ti)
IMPLICIT NONE

INTEGER(dp_int), INTENT(IN) :: nsti
REAL(dp_real), dimension(ndim4, nsti), INTENT(IN) :: b11_atms, o11_atms
REAL(dp_real),INTENT(IN) :: latt_vec(3,3), vec_unit(3)

INTEGER(dp_int), INTENT(INOUT) :: angle1_ti_o_ti(3,nsti)

REAL(dp_real) :: r1(3), r2(3), v1(3), v2(3), vo(3),  r1_norm, r2_norm

REAL(dp_real):: lz, dmin, ro1(3), ro2(3), r_dot, lsup(3)

INTEGER(dp_int) :: i, j,  k, no, ncount
INTEGER(dp_int) :: iz, l

! PRINT*, latt_vec

dmin=4.0/2.0 + 1.0
ncount = 0

OPEN(UNIT=71, FILE='angle_info_zaxis.dat', ACTION='WRITE')
WRITE(71,*) "***** O atoms along axis *****"
DO i = 1, nsti
 no = angle1_ti_o_ti(2, i)
 IF (no ==-1) THEN
    PRINT*, '==Program will give wrong results, All Ti do not have O atoms=='
    PRINT*, '**please check program**'
 ENDIF
 
 
! print*, no
 v1 = b11_atms(2:4, i )
 vo = o11_atms(2:4,no)

DO l =1, 3
 ro1(l) = v1(l) - vo(l)
ENDDO
r1_norm = SQRT(ro1(1)**2 + ro1(2)**2 + ro1(3)**2)
DO l =1, 3
  r1(l) = ro1(l)/r1_norm
ENDDO

DO j =1, nsti

 v2 = b11_atms(2:4, j)
 ro2(1)  = v2(1) -vo(1)
 ro2(2)  = v2(2) -vo(2)
 ro2(3)  = v2(3) -vo(3)

 r2_norm = SQRT(ro2(1)**2 + ro2(2)**2 + ro2(3)**2)
 r2(1) = ro2(1)/r2_norm
 r2(2) = ro2(2)/r2_norm
 r2(3) = ro2(3)/r2_norm

 r_dot = DOT_PRODUCT(r1,r2)
 IF (r_dot+1 <= 0.001 .AND. (r2_norm-dmin) <0.0) THEN 
     angle1_ti_o_ti(3, i) =j
    WRITE(71,*) i, no, j, r_dot
    ncount = ncount + 1
    IF (INT(r_dot) /=-1) PRINT*, "----WRONG---", r_dot 
 ENDIF


ENDDO
ENDDO


DO i = 1, nsti
k = angle1_ti_o_ti(3, i)
IF (k /=-1) cycle ! 

no = angle1_ti_o_ti(2, i)
v1 = b11_atms(2:4, i)
vo = o11_atms(2:4, no)

   
   DO iz = -1, 1, 2
     lz = DBLE(iz)
     lsup = (/latt_vec(1,1)*lz*vec_unit(1), latt_vec(2,2)*lz*vec_unit(2), latt_vec(3,3)*lz*vec_unit(3)/)

      ro1(1) = v1(1) -  (vo(1) + lsup(1))
      ro1(2) = v1(2) -  (vo(2) + lsup(2))
      ro1(3)  = v1(3) - (vo(3) + lsup(3))
       
     r1_norm = SQRT(ro1(1)**2 + ro1(2)**2 + ro1(3)**2)
     IF ((r1_norm-dmin) <0) THEN
            r1(1) = ro1(1)/r1_norm
            r1(2) = ro1(2)/r1_norm
            r1(3) = ro1(3)/r1_norm

            vo(1) = vo(1) + lsup(1)
            vo(2) = vo(2) + lsup(2)
            vo(3) = vo(3) +  lsup(3)
     ENDIF 
   ENDDO

  DO j =1, nsti
      IF(j==i) cycle
      v2 = b11_atms(2:4, j)
            
   DO iz = -1, 1, 2
      lz = DBLE(iz)
      lsup = (/latt_vec(1,1)*lz*vec_unit(1), latt_vec(2,2)*lz*vec_unit(2), latt_vec(3,3)*lz*vec_unit(3)/)
      ro2(1) =  v2(1) + lsup(1) - vo(1)
      ro2(2) =  v2(2) + lsup(2) - vo(2)
      ro2(3)  = v2(3) + lsup(3) - vo(3)
     r2_norm = SQRT(ro2(1)**2 + ro2(2)**2 + ro2(3)**2)
    IF ((r2_norm-dmin) <0) THEN
            r2(1) = ro2(1)/r2_norm
            r2(2) = ro2(2)/r2_norm
            r2(3) = ro2(3)/r2_norm 
            !print*, r2_norm
        angle1_ti_o_ti(3, i) =j
        WRITE(71,*) i, no, j, DOT_PRODUCT(r1,r2)
        ncount = ncount + 1
        IF (INT(DOT_PRODUCT(r1,r2)) /=-1) PRINT*, "----WRONG---", r_dot

    ENDIF

    ENDDO
 ENDDO

ENDDO
IF (ncount /= nsti) PRINT*, '-----WRONG number of angle along z-axis----'
!CLOSE(71)
END SUBROUTINE dm_Zint_angle

END MODULE data_cubic40
