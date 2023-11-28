
SUBROUTINE interaction_3atms_xy2(latt_sup, nsti, b_atms, o3_atms, angle_ti_o_ti)
   USE dtype, ONLY : dp_real, dp_int
   IMPLICIT NONE

   INTEGER(dp_int), INTENT(IN) :: nsti

   REAL(dp_real), INTENT(IN) ::  b_atms(nsti, 4), o3_atms(nsti,4), latt_sup(3,3)

   INTEGER(dp_int), INTENT(OUT) :: angle_ti_o_ti(nsti, 3)

   INTEGER(dp_int) ::  neigh_ab(nsti, 8), NeighAlongAxis(nsti)
   
   REAL (dp_real) :: dmin

   INTEGER(dp_int) :: l,k 


   DO l =1, nsti
    angle_ti_o_ti(l, 1) = l
      DO k =2, 3
      angle_ti_o_ti(l, k) = -1
      ENDDO 
   ENDDO
! Fe -O1
        dmin=4.0/2.0 + 1.0
        CALL neighbour_atm1_atm2(dmin, nsti, b_atms, o3_atms, neigh_ab)
        
!          DO l =1, nsti 
!                  PRINT*, 'Ti->O', l, (neigh_ab(l, k), k=1,6)
!          ENDDO
!          PRINT*,"----------"
     
       CALL projection_car_xy2(nsti, b_atms, o3_atms, neigh_ab, NeighAlongAxis)
       
        
 !  DO l =1, nsti 
 !     PRINT*,'3. Ti->O (+axis)', l,  (NeighAlongAxis(l))
 !  ENDDO
   PRINT*,"----------"
        
        DO l =1, nsti 
            if (NeighAlongAxis(l) /=0) angle_ti_o_ti(l,2) = NeighAlongAxis(l)
        ENDDO 
        
        CALL dm_int_angle(nsti, latt_sup, b_atms, o3_atms, angle_ti_o_ti)


END SUBROUTINE interaction_3atms_xy2
