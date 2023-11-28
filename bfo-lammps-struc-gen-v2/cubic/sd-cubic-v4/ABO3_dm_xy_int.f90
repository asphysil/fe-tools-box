SUBROUTINE dm_int_angle(nsti, latt_vec, b_atms, o_atms, angle_ti_o_ti)

        USE dtype, ONLY : dp_real, dp_int
  IMPLICIT NONE

  INTEGER(dp_int), INTENT(IN) :: nsti
  INTEGER(dp_int), INTENT(INOUT) :: angle_ti_o_ti(nsti, 3)
  REAL(dp_real), INTENT(IN) :: b_atms(nsti, 4), o_atms(nsti,4), latt_vec(3,3)


  REAL(dp_real) :: r1(3), r2(3), v1(3), v2(3), vo(3), r1_norm, r2_norm

REAL(dp_real):: lx, ly, ro1(3), ro2(3), r_dot, dmin

  INTEGER(dp_int) ::  i, j,  no, k
  INTEGER :: ix, iy, ncount

 ! PRINT*, latt_vec

 dmin=4.0/2.0 + 1.0
ncount = 0

 OPEN(UNIT=72, FILE='angle_info_xyaxis.dat', ACTION='WRITE')

 WRITE(72,*)"xy plane atoms Angle info"

  DO i = 1, nsti
     no = angle_ti_o_ti(i, 2)
     IF (no ==-1) THEN
   PRINT*, '==Program will give wrong results, All Ti do not have O atoms=='
   PRINT*, '**please check program**'
     ENDIF
     
     
    ! print*, no
     v1 = b_atms(i, 2:4)
     vo = o_atms(no, 2:4)

     ro1(1) = v1(1) -vo(1)
     ro1(2) = v1(2) -vo(2)
     ro1(3) = v1(3) -vo(3)
     r1_norm = SQRT(ro1(1)**2 + ro1(2)**2 + ro1(3)**2)

     r1(1) = ro1(1)/r1_norm
     r1(2) = ro1(2)/r1_norm
     r1(3) = ro1(3)/r1_norm

  DO j =1, nsti

     v2 = b_atms(j, 2:4)
     ro2(1)  = v2(1) - vo(1)
     ro2(2)  = v2(2) - vo(2)
     ro2(3)  = v2(3) - vo(3)
     r2_norm = SQRT(ro2(1)**2 + ro2(2)**2 + ro2(3)**2)
     
     r2(1) = ro2(1)/r2_norm
     r2(2) = ro2(2)/r2_norm
     r2(3) = ro2(3)/r2_norm

     r_dot = DOT_PRODUCT(r1,r2)

    !PRINT*, r_dot
    IF (r_dot+1 <= 0.001 .AND. (r2_norm-dmin) <0.0) THEN
     angle_ti_o_ti(i, 3) =j
      WRITE(72,*) i, no, j, r_dot
     ncount = ncount + 1
     IF (INT(r_dot) /=-1) PRINT*, "----WRONG---", r_dot
    ENDIF 

    ENDDO
ENDDO


DO i = 1, nsti
   k = angle_ti_o_ti(i, 3)
   IF (k /=-1) cycle
   
   no = angle_ti_o_ti(i, 2)
   v1 = b_atms(i, 2:4)
   vo = o_atms(no, 2:4)
      
      
       DO ix = -1, 1
          lx = DBLE(ix)
          DO iy = -1, 1 
              ly = DBLE(iy)

             ro1(1) = v1(1) -  (vo(1) + latt_vec(1,1)*lx)
             ro1(2) = v1(2) -  (vo(2) + latt_vec(2,2)*ly)
             ro1(3)  = v1(3) - vo(3) 

         r1_norm = SQRT(ro1(1)**2 + ro1(2)**2 + ro1(3)**2)
         IF ((r1_norm-dmin) <0) THEN
                r1(1) = ro1(1)/r1_norm
                r1(2) = ro1(2)/r1_norm
                r1(3) = ro1(3)/r1_norm
                !PRINT*, ncount, r1 ! DOT_PRODUCT(v1, vo)
                vo(1) = vo(1) + latt_vec(1,1)*lx
                vo(2) = vo(2) + latt_vec(2,2)*ly
                
         ENDIF
            ENDDO
       ENDDO

      DO j =1, nsti
          IF(j==i) cycle
          v2 = b_atms(j, 2:4)
       
      
       DO ix = -1, 1
          lx = DBLE(ix)
          DO iy = -1, 1
             ly = DBLE(iy)

          ro2(1) =  v2(1) + latt_vec(1,1)*lx -vo(1)
          ro2(2) =  v2(2) + latt_vec(2,2)*ly -vo(2)
          ro2(3)  = v2(3) - vo(3) 

         r2_norm = SQRT(ro2(1)**2 + ro2(2)**2 + ro2(3)**2)

        IF ((r2_norm-dmin) <0) THEN
                r2(1) = ro2(1)/r2_norm
                r2(2) = ro2(2)/r2_norm
                r2(3) = ro2(3)/r2_norm
            angle_ti_o_ti(i, 3) =j
            !ncount = ncount + 1
            !PRINT*, DOT_PRODUCT(r1,r2) 
            WRITE(72,*) i, no, j, DOT_PRODUCT(r1,r2)
            ncount = ncount + 1
           IF (INT(DOT_PRODUCT(r1,r2)) /=-1) PRINT*, "----WRONG---", r_dot
        ENDIF

        ENDDO
      ENDDO

     ENDDO

ENDDO
IF (ncount /= nsti) PRINT*, '-----WRONG number of angle along xy-axis----', ncount

END SUBROUTINE dm_int_angle
