
SUBROUTINE interaction_3atms_xy1(latt_sup, nsti, b_atms, o2_atms, angle_ti_o_ti)
   USE dtype, ONLY : dp_real, dp_int
   IMPLICIT NONE

   INTEGER(dp_int), INTENT(IN) :: nsti

   REAL(dp_real), INTENT(IN) ::  b_atms(nsti, 4), o2_atms(nsti,4), latt_sup(3,3)

   INTEGER(dp_int), INTENT(OUT) :: angle_ti_o_ti(nsti, 3)

   INTEGER(dp_int) ::  neigh_ab(nsti, 8), NeighAlongAxis(nsti), atmsij_tio(nsti)
   
   REAL (dp_real) :: dmin, vec_unit(3)

   INTEGER(dp_int) :: l, k,  nij_tio

   DO l =1, nsti
    angle_ti_o_ti(l, 1) = l
      DO k =2, 3
      angle_ti_o_ti(l, k) = -1
      ENDDO 
   ENDDO

! Fe -O1
        dmin=4.0/2.0 + 1.0
        CALL neighbour_atm1_atm2(dmin, nsti, b_atms, o2_atms, neigh_ab)
        
!          DO l =1, nsti 
!                  PRINT*, 'Ti->O', l, (neigh_ab(l, k), k=1,6)
!          ENDDO
!          PRINT*,"----------"
     
       CALL projection_car_xy1(nsti, b_atms, o2_atms, neigh_ab, NeighAlongAxis)
       
   !PRINT*,   'Ti->O did not find O-atom if it is 0 '   
   !DO l =1, nsti 
    !  PRINT*,'2. Ti->O (+axis) ', l,  (NeighAlongAxis(l))
   !ENDDO
   !PRINT*,"----------"
      
   ! Sorting which Ti did not find O-atoms
       nij_tio = 0
        DO l =1, nsti 
            if (NeighAlongAxis(l) /=0) THEN
                angle_ti_o_ti(l,2) = NeighAlongAxis(l)
                
            ELSE
                nij_tio = nij_tio + 1
                atmsij_tio(nij_tio) = l
            ENDIF
        ENDDO 

        ! DO l =1, nij_tio
       !      PRINT*, 'Ti did not find O -atoms', atmsij_tio(l)
        !  ENDDO

           !DO k =1, nsti
            !PRINT*, (angle_ti_o_ti(k,l), l=1,3)
           !ENDDO

        
        !PRINT*, '****Ti-----O*****'
        dmin = 4.0/2.0  + 1.0
        vec_unit = (/1.0, 0.0, 0.0/)
        CALL neighbour_axis(vec_unit, dmin, latt_sup, nsti, b_atms, o2_atms, nij_tio, atmsij_tio, NeighAlongAxis)
        
        DO l=1, nij_tio 
           k = atmsij_tio(l)
           angle_ti_o_ti(k,2) = NeighAlongAxis(k)
           !PRINT*, 'Ti->O (+axis)', k,  (NeighAlongAxis(k))
        ENDDO
        !PRINT*,"***----------***"
        !DO l =1, nsti 
        !    PRINT*,'2. Ti->O (+axis) ', l,  angle_ti_o_ti(l,2)
        ! ENDDO
        ! PRINT*,"***----------***"

        DO l =1, nsti
            IF (angle_ti_o_ti(k,2) ==-1) THEN 
                PRINT*, l, 'Ti atom did find O atom, Please check program'
            ENDIF

        ENDDO
        CALL dm_int_angle(nsti, latt_sup, b_atms, o2_atms, angle_ti_o_ti)
    
END SUBROUTINE interaction_3atms_xy1
