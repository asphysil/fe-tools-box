MODULE nearest_neighbour_list
use constant, only : dp 
use data_structure, only : ndim3, dp, latt_vec
private 
public :: nearest_neighbour_ab, nearest_neighbour_aa, &
nearest_neighbour_ab1, nearest_neighbour_aa1,&
nearest_neighbour_ab2, nearest_neighbour_aa2
contains 

SUBROUTINE nearest_neighbour_ab(n1st, n1, n2, dab, atm1, atm2, ab_loc)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n1st, n1, n2
  REAL(dp), INTENT(IN) :: dab, atm1(ndim3,n1), atm2(ndim3,n2)
  INTEGER, INTENT(OUT) :: ab_loc(n1st,n1)
  !
  INTEGER :: atm_info(27*n2)
  REAL(dp) ::  d(27*n2), dist, d_sq, dmax

  REAL(dp) :: v1(ndim3), vi(ndim3), vf(ndim3)
!
  INTEGER :: ncount 
  INTEGER :: i, j,  n, l1, l2, l3, loc, k

  dmax=100 ! angs

  DO i=1, n1 
      v1(1)=atm1(1, i)
      v1(2)=atm1(2, i)
      v1(3)=atm1(3, i)
      n=0
      DO j=1, n2
         vi(1)=atm2(1,j)
         vi(2)=atm2(2,j)
         vi(3)=atm2(3,j)

         DO l1=-1, 1
            DO l2 =-1, 1
               DO l3=-1, 1
                 n=n+1
                  vf(1)=vi(1)+ l1*latt_vec(1,1) !
                  vf(2)=vi(2)+ l2*latt_vec(2,2) !
                  vf(3)=vi(3)+ l3*latt_vec(3,3)

                  atm_info(n)=j
                  !atm_info(n,2)=vf(1)
                  !atm_info(n,3)=vf(2)
                  !atm_info(n,4)=vf(3)
                  d(n)=SQRT((v1(1)-vf(1))**2 + &
                            (v1(2)-vf(2))**2 + &
                            (v1(3)-vf(3))**2)
                 ENDDO
               ENDDO
             ENDDO
         ENDDO

          ! loc = MINLOC(d, DIM=1)
          ! d(loc) = 1000.0
!
           DO k=1, n1st
              loc=MINLOC(d, DIM=1)
               ab_loc(k, i)=INT(atm_info(loc))
               d(loc)=dmax ! to find new min value location
               ENDDO
    ENDDO 
  

END SUBROUTINE nearest_neighbour_ab

SUBROUTINE nearest_neighbour_aa(n1st, n1, n2, dab, atm1, atm2, ab_loc)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n1st, n1, n2
  REAL(dp), INTENT(IN) :: dab, atm1(ndim3,n1), atm2(ndim3,n2)
  INTEGER, INTENT(OUT) :: ab_loc(n1st,n1)
  !
  INTEGER :: atm_info(27*n2)
  REAL(dp) ::  d(27*n2), dist, d_sq, dmax

  REAL(dp) :: v1(ndim3), vi(ndim3), vf(ndim3)
!
  INTEGER :: ncount 
  INTEGER :: i, j,  n, l1, l2, l3, loc, k

  dmax=100 ! angs

  DO i=1, n1 
      v1(1)=atm1(1, i)
      v1(2)=atm1(2, i)
      v1(3)=atm1(3, i)
      n=0
      DO j=1, n2
         vi(1)=atm2(1,j)
         vi(2)=atm2(2,j)
         vi(3)=atm2(3,j)

         DO l1=-1, 1
            DO l2 =-1, 1
               DO l3=-1, 1
                 n=n+1
                  vf(1)=vi(1)+ l1*latt_vec(1,1) !
                  vf(2)=vi(2)+ l2*latt_vec(2,2) !
                  vf(3)=vi(3)+ l3*latt_vec(3,3)

                  atm_info(n)=j
                  !atm_info(n,2)=vf(1)
                  !atm_info(n,3)=vf(2)
                  !atm_info(n,4)=vf(3)
                  d(n)=SQRT((v1(1)-vf(1))**2 + &
                            (v1(2)-vf(2))**2 + &
                            (v1(3)-vf(3))**2)
                 ENDDO
               ENDDO
             ENDDO
         ENDDO

           loc = MINLOC(d, DIM=1)
           d(loc) = 1000.0
!
           DO k=1, n1st
              loc=MINLOC(d, DIM=1)
               ab_loc(k, i)=INT(atm_info(loc))
               d(loc)=dmax ! to find new min value location
               ENDDO
    ENDDO 
  

END SUBROUTINE nearest_neighbour_aa

SUBROUTINE nearest_neighbour_ab1(n1st, n1, n2, dab, atm1, atm2, ab_loc)
 ! dab is minimum a-b distance
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n1st, n1, n2
    REAL(dp), INTENT(IN) :: dab, atm1(ndim3,n1), atm2(ndim3,n2)
    INTEGER, INTENT(OUT) :: ab_loc(n1st,n1)
    !
    INTEGER :: atm_info(27)
    REAL(dp) ::  d(27), dist, d_sq, dmax
  
    REAL(dp) :: v1(ndim3), vi(ndim3), vf(ndim3)
  !
    INTEGER :: ncount, nn  
    INTEGER :: i, j,  n, l1, l2, l3, loc
  
     dmax=100 ! angs
     
  d_sq = dab*dab  
  dmax=100 ! angs

  DO i=1, n1 
      v1(1)=atm1(1, i)
      v1(2)=atm1(2, i)
      v1(3)=atm1(3, i)
      
      ncount = 0 
      DO j=1, n2
         vi(1)=atm2(1,j)
         vi(2)=atm2(2,j)
         vi(3)=atm2(3,j)

        if (ncount == n1st) exit

        dist = (v1(1)-vi(1))**2 + &
                (v1(2)-vi(2))**2 + &
                (v1(3)-vi(3))**2

        if (dist < d_sq) then
            ncount = ncount + 1
            ab_loc(ncount, i)= j 
        cycle
        else
          nn = 0 
          DO l1=-1, 1
            DO l2 =-1, 1
               DO l3=-1, 1
                 nn=nn+1
                  vf(1)=vi(1)+ l1*latt_vec(1,1) !
                  vf(2)=vi(2)+ l2*latt_vec(2,2) !
                  vf(3)=vi(3)+ l3*latt_vec(3,3)

                  atm_info(nn)=j
                  !atm_info(n,2)=vf(1)
                  !atm_info(n,3)=vf(2)
                  !atm_info(n,4)=vf(3)
                  d(nn)=(v1(1)-vf(1))**2 + &
                            (v1(2)-vf(2))**2 + &
                            (v1(3)-vf(3))**2
                 ENDDO
               ENDDO
             ENDDO
        endif 
        loc=MINLOC(d, DIM=1) 
        dist = d(loc) 
        if (dist < d_sq) then
          ncount = ncount + 1
          ab_loc(ncount, i)= INT(atm_info(loc))
        endif 
      ENDDO
    ENDDO 
  
  END SUBROUTINE nearest_neighbour_ab1


SUBROUTINE nearest_neighbour_ab2(n1st, n1, n2, dab, atm1, atm2, ab_loc)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n1st, n1, n2
  REAL(dp), INTENT(IN) :: dab, atm1(ndim3,n1), atm2(ndim3,n2)
  INTEGER, INTENT(OUT) :: ab_loc(n1st,n1)
  !
  INTEGER :: atm_info(27), index_loc(500)
  REAL(dp) ::  d(27), d_temp(500), dist, d_sq, dmax

  REAL(dp) :: v1(ndim3), vi(ndim3), vf(ndim3)
!
  INTEGER :: ncount
  INTEGER :: ii, j,  nn, l1, l2, l3, loc, k

   dmax=100 ! angs

d_sq = dab*dab

   DO ii=1, n1

       v1(1)=atm1(1,ii)
       v1(2)=atm1(2,ii)
       v1(3)=atm1(3,ii)
       ncount =0
      

       DO j=1, n2
        !print*, ii, j 

          vi(1)=atm2(1,j)
          vi(2)=atm2(2,j)
          vi(3)=atm2(3,j)

          dist = (v1(1)-vi(1))**2 + &
                 (v1(2)-vi(2))**2 + &
                 (v1(3)-vi(3))**2

          if (dist < d_sq) then
            ncount = ncount + 1
            d_temp(ncount) = dist 
            index_loc(ncount) = j !
          ELSE
            nn=0
            DO l1=-1, 1
              DO l2 =-1, 1
                 DO l3=-1, 1
                   nn=nn+1
                    vf(1)=vi(1)+ DBLE(l1)*latt_vec(1,1) !
                    vf(2)=vi(2)+ DBLE(l2)*latt_vec(2,2) !
                    vf(3)=vi(3)+ DBLE(l3)*latt_vec(3,3)
                    !print*, vf(1), vf(2), vf(3) 
          !          !atm_info(n)=j
                    d(nn)=(v1(1)-vf(1))*(v1(1)-vf(1)) + &
                          (v1(2)-vf(2))*(v1(2)-vf(2)) + &
                          (v1(3)-vf(3))*(v1(3)-vf(3))
                   ENDDO
                 ENDDO
            ENDDO
            loc=MINLOC(d, DIM=1)
            dist = d(loc)
          endif

          if (dist < d_sq) then
              ncount = ncount + 1
              d_temp(ncount) = dist 
              index_loc(ncount) = j 
          endif

      ENDDO

  DO k=1, n1st
    loc = MINLOC(d_temp(1:ncount), DIM=1)
    ab_loc(k,ii) = index_loc(loc)
    d_temp(loc) = 1000.0
  ENDDO
  !print*,ii, n1, ncount !,  m, d_sq
  ENDDO

END SUBROUTINE nearest_neighbour_ab2

SUBROUTINE nearest_neighbour_aa1(n1st, n1, n2, dab, atm1, atm2, ab_loc)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n1st, n1, n2
    REAL(dp), INTENT(IN) :: dab, atm1(ndim3,n1), atm2(ndim3,n2)
    INTEGER, INTENT(OUT) :: ab_loc(n1st,n1)
    !
    INTEGER :: atm_info(27)
    REAL(dp) ::  d(27), dist, d_sq, dmax

    REAL(dp) :: v1(ndim3), vi(ndim3), vf(ndim3)
  !
    INTEGER :: ncount, nn 
    INTEGER :: i, j,  n, l1, l2, l3, loc, k

     dmax=100 ! angs
 
  d_sq = dab*dab
  DO i=1, n1 
    v1(1)=atm1(1, i)
    v1(2)=atm1(2, i)
    v1(3)=atm1(3, i)
    n=0
    ncount = 0 
    DO j=1, n2
      if(j==i) cycle  
      !print*, j, d_sq, dist, n2, ncount  
       vi(1)=atm2(1,j)
       vi(2)=atm2(2,j)
       vi(3)=atm2(3,j)

      if (ncount == n1st) exit

      dist = (v1(1)-vi(1))**2 + &
              (v1(2)-vi(2))**2 + &
              (v1(3)-vi(3))**2

      if (dist < d_sq) then
          ncount = ncount + 1
          ab_loc(ncount, i)= j
          
      cycle
      else
        !print*, "ok"
        nn = 0 
        DO l1=-1, 1
          DO l2 =-1, 1
             DO l3=-1, 1
               nn=nn+1
                vf(1)=vi(1)+ l1*latt_vec(1,1) !
                vf(2)=vi(2)+ l2*latt_vec(2,2) !
                vf(3)=vi(3)+ l3*latt_vec(3,3)

                atm_info(nn)=j
                !atm_info(n,2)=vf(1)
                !atm_info(n,3)=vf(2)
                !atm_info(n,4)=vf(3)
                d(nn)=(v1(1)-vf(1))**2 + &
                          (v1(2)-vf(2))**2 + &
                          (v1(3)-vf(3))**2
               ENDDO
             ENDDO
           ENDDO
      endif 
      loc=MINLOC(d, DIM=1) 
      dist = d(loc) 
      if (dist < d_sq) then
        ncount = ncount + 1
        ab_loc(ncount, i)= INT(atm_info(loc))
      endif 
    ENDDO
  ENDDO 
  END SUBROUTINE nearest_neighbour_aa1


SUBROUTINE nearest_neighbour_aa2(n1st, n1, n2, dab, atm1, atm2, ab_loc)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n1st, n1, n2
  REAL(dp), INTENT(IN) :: dab, atm1(ndim3,n1), atm2(ndim3,n2)
  INTEGER, INTENT(OUT) :: ab_loc(n1st,n1)
  !
  INTEGER :: atm_info(27), index_loc(500)
  REAL(dp) ::  d(27), d_temp(500), dist, d_sq, dmax

  REAL(dp) :: v1(ndim3), vi(ndim3), vf(ndim3)
!
  INTEGER :: ncount, m
  INTEGER :: ii, j,  nn, l1, l2, l3, loc, k

   dmax=100 ! angs

d_sq = dab*dab
!print*, d_sq
!print*, latt_vec(1,1), latt_vec(2,2), latt_vec(3,3)

   DO ii=1, n1

       v1(1)=atm1(1,ii)
       v1(2)=atm1(2,ii)
       v1(3)=atm1(3,ii)
       ncount =0
       !print*, n1, n2, v1(1),v1(2),v1(3) 
       !m=0

       DO j=1, n2
        !print*, ii, j 
         if (j==ii) CYCLE
        ! m = m + 1 

          vi(1)=atm2(1,j)
          vi(2)=atm2(2,j)
          vi(3)=atm2(3,j)

          dist = (v1(1)-vi(1))**2 + &
                 (v1(2)-vi(2))**2 + &
                 (v1(3)-vi(3))**2

          if (dist < d_sq) then
            ncount = ncount + 1
            d_temp(ncount) = dist 
            index_loc(ncount) = j !
          ELSE
            nn=0
            DO l1=-1, 1
              DO l2 =-1, 1
                 DO l3=-1, 1
                   nn=nn+1
                    vf(1)=vi(1)+ DBLE(l1)*latt_vec(1,1) !
                    vf(2)=vi(2)+ DBLE(l2)*latt_vec(2,2) !
                    vf(3)=vi(3)+ DBLE(l3)*latt_vec(3,3)
                    !print*, vf(1), vf(2), vf(3) 
          !          !atm_info(n)=j
                    d(nn)=(v1(1)-vf(1))*(v1(1)-vf(1)) + &
                          (v1(2)-vf(2))*(v1(2)-vf(2)) + &
                          (v1(3)-vf(3))*(v1(3)-vf(3))
                   ENDDO
                 ENDDO
            ENDDO
            loc=MINLOC(d, DIM=1)
            dist = d(loc)
          endif

          if (dist < d_sq) then
              ncount = ncount + 1
              d_temp(ncount) = dist 
              index_loc(ncount) = j 
          endif

      ENDDO

  DO k=1, n1st
    loc = MINLOC(d_temp(1:ncount), DIM=1)
    ab_loc(k,ii) = index_loc(loc)
    d_temp(loc) = 1000.0
  ENDDO
  !print*,ii, n1, ncount !,  m, d_sq
  ENDDO

END SUBROUTINE nearest_neighbour_aa2




SUBROUTINE nearest_neighbour_fe1fe6(n1st, nfe, dab, atm1,ab_loc, fe1fe6_coord)
! n1st = number of nearest neighbour atoms
! nfe = number of Fe atoms
! dab = minimum distance between Fe1-Fe2 
! atm1 = Fe atoms coordinates 
! ab_loc = index of  nearest neighbour six atoms 
! fe1fe6_coord = output of six atoms coordinates

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n1st, nfe
  REAL(dp), INTENT(IN) :: dab, atm1(ndim3,nfe)
  INTEGER, INTENT(IN) :: ab_loc(n1st,nfe)
  REAL(dp), INTENT(OUT) :: fe1fe6_coord(3,n1st,nfe)
  !
  REAL(dp) :: atm_info(3, 27), temp_coord(ndim3, n1st, nfe)
  REAL(dp) ::  d(27), dist, d_sq, dmax

  REAL(dp) :: v1(ndim3), vi(ndim3), vf(ndim3)
!
  INTEGER :: ncount, m
  INTEGER :: ii, j,  nn, l1, l2, l3, loc, k

   dmax=100 ! angs

d_sq = dab*dab
!print*, d_sq
!print*, latt_vec(1,1), latt_vec(2,2), latt_vec(3,3)
DO ii =1, nfe 
  DO j = 1, n1st
    m = ab_loc(j,ii)
    temp_coord(1:3,j,ii) = atm1(1:3,m)
  ENDDO
ENDDO

   DO ii=1, nfe
       v1(1:3)=atm1(1:3,ii)
       DO j=1, n1st 
          vi(1:3)=temp_coord(1:3, j, ii)

          dist = (v1(1)-vi(1))**2 + &
                 (v1(2)-vi(2))**2 + &
                 (v1(3)-vi(3))**2

          if (dist < d_sq) then
            fe1fe6_coord(1:3,j,ii) =vi(1:3)
          ELSE
            nn=0
            DO l1=-1, 1
              DO l2 =-1, 1
                 DO l3=-1, 1
                   nn=nn+1
                    vf(1)=vi(1)+ DBLE(l1)*latt_vec(1,1) !
                    vf(2)=vi(2)+ DBLE(l2)*latt_vec(2,2) !
                    vf(3)=vi(3)+ DBLE(l3)*latt_vec(3,3)
                    !print*, vf(1), vf(2), vf(3)

                   atm_info(1,nn)=vf(1)
                   atm_info(2,nn)=vf(2)
                   atm_info(3,nn)=vf(3)

                    d(nn)=(v1(1)-vf(1))*(v1(1)-vf(1)) + &
                          (v1(2)-vf(2))*(v1(2)-vf(2)) + &
                          (v1(3)-vf(3))*(v1(3)-vf(3))
                   ENDDO
                 ENDDO
            ENDDO
            loc=MINLOC(d, DIM=1)
           fe1fe6_coord(1:3, j, ii) = atm_info(1:3, loc)
          endif

      ENDDO
  ENDDO

END SUBROUTINE nearest_neighbour_fe1fe6


SUBROUTINE angle_coord_o1feo2(no, nfe, dab, atm_O, atm_fe, angle_index, coord_fe1ofe2)
! Number of angle = number of O atoms (no)
!  no = number of o atoms
! dab = minimum distance between O-Fe
! atm_O = O atoms coordinates 
! atm_fe = Fe atoms coordinates 
! angle_index = index of three atoms 
! coord_fe1Ofe2 = output of three atoms coordinates
!  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: no, nfe
  REAL(dp), INTENT(IN) :: dab, atm_O(ndim3,no), atm_fe(ndim3, nfe)
  INTEGER, INTENT(IN) :: angle_index(3, no)
  REAL(dp), INTENT(OUT) :: coord_fe1ofe2(ndim3,3,no)
  !
  REAL(dp) :: atm_info(3, 27), temp_coord(3, 3, no)
  REAL(dp) ::  d(27), dist1, dist2, d_sq, dmax

  REAL(dp) :: v1(ndim3), w(ndim3), v2(ndim3), vf(ndim3)
!
  INTEGER ::  m1, m2, m3
  INTEGER :: ii, j,  nn, l1, l2, l3, loc

   dmax=100 ! angs

d_sq = dab*dab
!print*, d_sq
!print*, latt_vec(1,1), latt_vec(2,2), latt_vec(3,3)

DO ii =1, no 
    m1 = angle_index(1,ii) ! Fe1
    m2 = angle_index(2,ii) ! O
    m3 = angle_index(3,ii) ! Fe2

    temp_coord(1:3, 1, ii) = atm_fe(1:3,m1)
    temp_coord(1:3, 2, ii) = atm_O(1:3,m2)
    temp_coord(1:3, 3, ii) = atm_fe(1:3,m3)
  ENDDO

   DO ii=1, no

       v1(1:3)=temp_coord(1:3, 1, ii) ! Fe1
       w(1:3)=temp_coord(1:3,  2, ii) ! O
       v2(1:3)=temp_coord(1:3, 3, ii) ! Fe2

        dist1 = (w(1)-v1(1))**2 + &
               (w(2)-v1(2))**2 + &
               (w(3)-v1(3))**2

        dist2 = (w(1)-v2(1))**2 + &
                (w(2)-v2(2))**2 + &
                (w(3)-v2(3))**2

        coord_fe1ofe2(1:3,2,ii) = w(1:3) ! O

! minimum distance (O1-Fe1) coordinates of Fe1
          if (dist1 < d_sq) then
          coord_fe1ofe2(1:3,1,ii) =v1(1:3) ! Fe1
          ELSE
            nn=0
            DO l1=-1, 1
              DO l2 =-1, 1
                 DO l3=-1, 1
                   nn=nn+1
                    vf(1)=v1(1)+ DBLE(l1)*latt_vec(1,1) !
                    vf(2)=v1(2)+ DBLE(l2)*latt_vec(2,2) !
                    vf(3)=v1(3)+ DBLE(l3)*latt_vec(3,3)
                    !print*, vf(1), vf(2), vf(3)

                   atm_info(1,nn)=vf(1)
                   atm_info(2,nn)=vf(2)
                   atm_info(3,nn)=vf(3)

                    d(nn)=(w(1)-vf(1))**2 + &
                          (w(2)-vf(2))**2 + &
                          (w(3)-vf(3))**2
                   ENDDO
                 ENDDO
            ENDDO
            loc=MINLOC(d, DIM=1)
          coord_fe1ofe2(1:3,1,ii) =atm_info(1:3, loc) ! Fe1
          endif
! minimum distance (O1-Fe2) coordinates of Fe2 
      if (dist2 < d_sq) then
        coord_fe1ofe2(1:3,3,ii) = v2(1:3) ! Fe1
          ELSE
            nn=0
            DO l1=-1, 1
              DO l2 =-1, 1
                 DO l3=-1, 1
                   nn=nn+1
                    vf(1)=v2(1)+ DBLE(l1)*latt_vec(1,1) !
                    vf(2)=v2(2)+ DBLE(l2)*latt_vec(2,2) !
                    vf(3)=v2(3)+ DBLE(l3)*latt_vec(3,3)
                    !print*, vf(1), vf(2), vf(3)

                   atm_info(1,nn)=vf(1)
                   atm_info(2,nn)=vf(2)
                   atm_info(3,nn)=vf(3)

                    d(nn)=(w(1)-vf(1))**2 + &
                          (w(2)-vf(2))**2 + &
                          (w(3)-vf(3))**2
                   ENDDO
                 ENDDO
            ENDDO
            loc=MINLOC(d, DIM=1)
          coord_fe1ofe2(1:3,3,ii) =atm_info(1:3, loc) ! Fe1
          endif
  ENDDO

END SUBROUTINE angle_coord_o1feo2

  END MODULE nearest_neighbour_list
  
