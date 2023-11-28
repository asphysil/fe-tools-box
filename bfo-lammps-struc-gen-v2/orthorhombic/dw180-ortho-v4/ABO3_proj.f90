SUBROUTINE projection_car_axis(vec_unit, nsti, atm1, atm2, neigh_ab, NeighAlongAxis)

        USE dtype, ONLY : dp_real, dp_int
  IMPLICIT NONE 

  INTEGER(dp_int), INTENT(IN) :: nsti, neigh_ab(nsti, 8)

        REAL(dp_real), INTENT(IN) :: atm1(nsti, 4), atm2(nsti,4), vec_unit(3)

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
    v1 = atm1(i, 2:4)
    ncout = 0
    DO j = 1, 8
      m = neigh_ab(i,j)
      IF(m==-1) CYCLE
        v2 = atm2(m, 2:4)

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
