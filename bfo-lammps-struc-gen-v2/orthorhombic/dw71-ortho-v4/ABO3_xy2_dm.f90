
SUBROUTINE interaction_3atms_xy2(latt_sup, nsti, b1_atms, o33_atms, angle1_ti_o_ti)

   USE my_constants, ONLY : dp_real, dp_int
   USE data_cubic40, ONLY : neighbour_atm1_atm2,neighbour_axis,&
                            projection_car_xy2, dm_Zint_angle, dm_int_angle

   IMPLICIT NONE

   INTEGER(dp_int), INTENT(IN) :: nsti

   REAL(dp_real), INTENT(IN) ::  b1_atms(4, nsti), o33_atms(4, nsti), latt_sup(3,3)

   INTEGER(dp_int), INTENT(OUT) :: angle1_ti_o_ti(3, nsti)

   INTEGER(dp_int) ::  neigh_ab(8, nsti), NeighAlongAxis(nsti)
   
   REAL (dp_real) :: dmin

   INTEGER(dp_int) :: l,k 


   DO l =1, nsti
    angle1_ti_o_ti(1, l) = l
      DO k =2, 3
      angle1_ti_o_ti(k, l) = -1
      ENDDO 
   ENDDO
! Fe -O1
        dmin=4.0/2.0 + 1.0
        CALL neighbour_atm1_atm2(dmin, nsti, b1_atms, o33_atms, neigh_ab)
        
!          DO l =1, nsti 
!                  PRINT*, 'Ti->O', l, (neigh_ab(l, k), k=1,6)
!          ENDDO
!          PRINT*,"----------"
     
       CALL projection_car_xy2(nsti, b1_atms, o33_atms, neigh_ab, NeighAlongAxis)
       
        
 !  DO l =1, nsti 
 !     PRINT*,'3. Ti->O (+axis)', l,  (NeighAlongAxis(l))
 !  ENDDO
   PRINT*,"----------"
        
        DO l =1, nsti 
            if (NeighAlongAxis(l) /=0) angle1_ti_o_ti(2,l) = NeighAlongAxis(l)
        ENDDO 
        
        CALL dm_int_angle(nsti, latt_sup, b1_atms, o33_atms, angle1_ti_o_ti)


END SUBROUTINE interaction_3atms_xy2
