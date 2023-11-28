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


  END MODULE nearest_neighbour_list
  
