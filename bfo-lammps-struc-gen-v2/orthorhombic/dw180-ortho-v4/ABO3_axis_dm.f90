
SUBROUTINE interaction_3atms(vec_unit, latt_sup, nsti, b1_atms, o11_atms, angle1_ti_o_ti)

   USE my_constants, ONLY : dp_real, dp_int
   USE data_cubic40, ONLY : neighbour_atm1_atm2,neighbour_axis,&
               projection_car_axis,dm_Zint_angle

   IMPLICIT NONE

   INTEGER(dp_int), INTENT(IN) :: nsti

   REAL(dp_real), INTENT(IN) :: vec_unit(3), b1_atms(4, nsti), o11_atms(4, nsti), latt_sup(3,3)

   INTEGER(dp_int), INTENT(OUT) :: angle1_ti_o_ti(3, nsti)

   INTEGER(dp_int) ::  neigh_ab(8, nsti), atmsij(nsti), NeighAlongAxis(nsti)
   
   REAL (dp_real) :: dmin

   INTEGER(dp_int) :: l, k, nij

   DO l =1, nsti
    angle1_ti_o_ti(1, l) = l
      DO k =2, 3
      angle1_ti_o_ti(k, l) = -1
      ENDDO 
   ENDDO
! Fe -O1
        dmin=4.0/2.0 + 1.0
        CALL neighbour_atm1_atm2(dmin, nsti, b1_atms, o11_atms, neigh_ab)
        
         !  DO l =1, nsti 
         !          PRINT*, 'Ti->O', l, (neigh_ab(l, k), k=1,6)
         !  ENDDO
         !  PRINT*,"----------"
        CALL projection_car_axis(vec_unit, nsti, b1_atms, o11_atms, neigh_ab, NeighAlongAxis)
        
   ! DO l =1, nsti 
   !    PRINT*,'Ti->O (+axis)', l,  (NeighAlongAxis(l))
   ! ENDDO
   ! PRINT*,"----------"
        nij =0 
        DO l =1, nsti 
            if (NeighAlongAxis(l) /=0) THEN
               angle1_ti_o_ti(2, l) = NeighAlongAxis(l)
            ELSE
              nij = nij + 1
              atmsij(nij) = l
            ENDIF 
        ENDDO 
        
      !   ! Fe1-Fe2
      !   dmin=4.0 + 1.0
      !   CALL neighbour_atm1_atm2(dmin, nsti, b_atms, b_atms, neigh_ab)
        
      ! !   DO l =1, nsti 
      ! !           PRINT*, 'Ti->Ti', l, (neigh_ab(l, k), k=1,6)
      ! !   ENDDO
      ! !   PRINT*,"----------"
      !   CALL projection_car_axis(vec_unit, nsti, b_atms, b_atms, neigh_ab, NeighAlongAxis)
        
      ! ! DO l =1, nsti 
      ! ! PRINT*, 'Ti->Ti (+axis)', l,  (NeighAlongAxis(l))
      ! ! ENDDO
      ! ! PRINT*,"----------"
        
      !   DO l =1, nsti 
      !       angle_ti_o_ti(l,1) = l 
      !       if (NeighAlongAxis(l) /=0) angle_ti_o_ti(l,3) = NeighAlongAxis(l)
      !   ENDDO 
        
        
        
        ! nij =0 
        ! DO l =1, nsti 
        !     IF (NeighAlongAxis(l) ==0) THEN
        !         nij = nij + 1
        !         atmsij(nij) = l 
        !     ENDIF
        ! ENDDO
        
        !  DO i =1, naij
        !      PRINT*, i, atmsij(i)
        !  ENDDO
        
        ! dmin = 4.0 +1.0
        ! CALL neighbour_axis(vec_unit, dmin, latt_sup, nsti, b_atms, b_atms, nij, atmsij, NeighAlongAxis)
      
        ! DO l =1, nij 
        ! k = atmsij(l)
        ! angle_ti_o_ti(k,3) = NeighAlongAxis(k)
        ! ENDDO 
        
        dmin = 4.0/2.0  +1.0
        CALL neighbour_axis(vec_unit, dmin, latt_sup, nsti, b1_atms, o11_atms, nij, atmsij, NeighAlongAxis)
        
        DO l=1, nij 
           k = atmsij(l)
           angle1_ti_o_ti(2, k) = NeighAlongAxis(k)
        ENDDO 

       ! DO l =1, nsti 
       !   PRINT*, '1. Ti->O (+axis)', l, angle_ti_o_ti(l,2)
       ! ENDDO
        

      CALL dm_Zint_angle(nsti, latt_sup, b1_atms, o11_atms, angle1_ti_o_ti)

END SUBROUTINE interaction_3atms
