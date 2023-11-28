SUBROUTINE fe_atoms_afe(nsti, atm1,  fe_mag )
         USE my_constants, ONLY : dp_real, dp_int
         USE data_cubic40, ONLY : neighbour_atm1_atm2
         IMPLICIT NONE
         INTEGER(dp_int), INTENT(IN) :: nsti 
         REAL(dp_real), INTENT(IN) :: atm1(4,nsti)
         INTEGER(dp_int), INTENT(OUT) :: fe_mag(nsti)
         
         INTEGER(dp_int) :: i, j, m, n, flag(nsti), neigh_ab(8,nsti)
         REAL(dp_real) :: dmin
         
         DO i =1, nsti 
            fe_mag(i)=-1
            flag(i)=1
         ENDDO 
         
         dmin = 4.0 + 1.0
         CALL neighbour_atm1_atm2(dmin, nsti, atm1, atm1, neigh_ab)

         DO i = 1, nsti 
            n = flag(i)
            IF (n==0) CYCLE  
         
            DO j = 1, 8 
           m = neigh_ab(j, i)
           IF (m==-1)CYCLE 
             
                 IF ( fe_mag(m) == -1) THEN 
                    IF (flag(m) == 1) THEN 
                   fe_mag(m) = 1
             flag(m) = 0
         
             ENDIF
              ENDIF
           ENDDO
         ENDDO

 j = 0
 DO i =1, nsti
     IF ( fe_mag(i) ==1) THEN
          j = j + 1
      ENDIF
 ENDDO
 PRINT*, '****Fe UP =', j, '*****Fe doun', nsti -j

END SUBROUTINE fe_atoms_afe
