
SUBROUTINE neighbour_atm1_atm2(dmin, nsti, atm1, atm2, neigh_ab)

        USE dtype, ONLY : dp_real, dp_int
    IMPLICIT NONE 
    INTEGER(dp_int), INTENT(IN) :: nsti
        REAL(dp_real), INTENT(IN) :: dmin, atm1(nsti, 4), atm2(nsti,4)

    INTEGER(dp_int), INTENT(OUT) :: neigh_ab(nsti,8)

    REAL(dp_real) :: v1(3), v2(3), dist

    INTEGER(dp_int) :: i, j

    INTEGER :: m


    DO i =1, nsti
        DO j = 1, 8
            neigh_ab(i,j) = -1
        ENDDO
    ENDDO

    DO i =1, nsti 
        v1= atm1(i,2:4)
        m =0

        DO j =1, nsti 
            v2 = atm2(j,2:4)

            dist = SQRT(( v1(1) - v2(1))**2 + &
                   (v1(2) - v2(2))**2   + &
         (v1(3) - v2(3))**2)

            IF (dist < 1.0) THEN
                CYCLE
            ENDIF  

            if ( (dist-dmin) < 0.0) THEN 
                !print*, dist
                m = m +1
                neigh_ab(i, m) = j
                ENDIF

        ENDDO 
    ENDDO

END SUBROUTINE neighbour_atm1_atm2


SUBROUTINE neighbour_axis(vec_unit, dmin, latt_sup, nsti, atm1, atm2, nij, atmsij, neighij)
    
        USE dtype, ONLY : dp_real, dp_int

        IMPLICIT NONE 
        INTEGER(dp_int), INTENT(IN) :: nsti, nij
        REAL(dp_real), INTENT(IN) :: dmin, atm1(nsti, 4), atm2(nsti,4), latt_sup(3,3), vec_unit(3)
        INTEGER(dp_int), INTENT(IN) :: atmsij(nsti)
    INTEGER(dp_int), INTENT(OUT) :: neighij(nsti)

    REAL(dp_real) :: v1(3), v2(3), w(3), dist, T(3), sl

    INTEGER(dp_int) :: i, j, m, k



    DO i =1, nsti
        neighij(i) = 0
    ENDDO

    DO i =1, nij
        m = atmsij(i)
        v1= atm1(m,2:4)

        DO j =1, nsti 
            v2 = atm2(j,2:4)

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
