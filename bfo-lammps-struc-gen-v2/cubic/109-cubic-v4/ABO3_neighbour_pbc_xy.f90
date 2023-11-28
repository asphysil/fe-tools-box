

SUBROUTINE neighbour_xy_pbc(dmin, latt_sup, nsti, atm1, atm2, atm3, nij, atmsij, neigh_tio, neighij)
    
        USE dtype, ONLY : dp_real, dp_int

        IMPLICIT NONE 
        INTEGER(dp_int), INTENT(IN) :: nsti, nij
        REAL(dp_real), INTENT(IN) :: dmin, atm1(nsti, 4), atm2(nsti,4), atm3(nsti,4), latt_sup(3,3)
        INTEGER(dp_int), INTENT(IN) :: atmsij(nsti),neigh_tio(nsti,3)
    INTEGER(dp_int), INTENT(OUT) :: neighij(nsti)

    REAL(dp_real) :: v1(3), v2(3), w(3), v_oxy(3), dist, dist_tio, T(2), slx, sly,  dmin_tio

    INTEGER(dp_int) :: i, j, m,  m_oxy
    INTEGER :: kx, ky, mm


    dmin_tio = dmin/2  

    DO i =1, nsti
        neighij(i) = 0
    ENDDO

    DO i =1, nij
        m = atmsij(i)
        v1= atm1(m,2:4)

        m_oxy = neigh_tio(m,2)
        !print*, m, m_oxy

        v_oxy = atm3(m_oxy,2:4)

        mm = 0

        DO j =1, nsti 
            v2 = atm2(j,2:4)

            DO kx = -1, 1
                DO ky = -1, 1
               
               IF ( (kx==0) .AND. (ky==0)) CYCLE 
                !print*,k
               slx = DBLE(kx)
               sly = DBLE(ky)

               T(1) =  latt_sup(1,1)*slx
               T(2) =  latt_sup(2,2)*sly
               

               w(1) = v2(1) + T(1)
               w(2) = v2(2) + T(2)
               w(3) = v2(3) 

               dist = SQRT((v1(1)-w(1))**2 + (v1(2)-w(2))**2 + (v1(3)-w(3))**2)
               !print*, dist
               
               if ( (dist-dmin) < 0.0) THEN 
                  !print*, m, j, dist, kx, ky        
                 dist_tio =  SQRT((v_oxy(1)-w(1))**2 + (v_oxy(2)-w(2))**2 + (v_oxy(3)-w(3))**2)
                 if (( dist_tio-dmin_tio)<0.0) THEN
                         !print*, dist, dist_tio
                         mm = mm + 1
                         neighij(m) = j
                !neigh_coord(i, mm, 1:3) = w(1:3) 
                    ENDIF
              ENDIF
           ENDDO
        ENDDO

               IF (mm>1) PRINT*, "mm value greater than 1"

        ENDDO 
    ENDDO

!    DO i = 1, nij
!       PRINT*, (neigh_index(i, j), j=1,3)
!    ENDDO

END SUBROUTINE neighbour_xy_pbc
