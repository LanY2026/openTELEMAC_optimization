!     Authored by Ziyi Huang, January 2026
!     Computation of infiltration water node by node, in coupling modeling of TELEMAC-2D (version tag v9p0r0) and modified Horton infiltration model. SMH is revised if there is infiltration water.

!     ******************************************************************
!     This subroutine should be called by a subroutine that sets SMH, which generally is PROSOU.
!     Execution of pre-computation is required before this subroutine is called.
!     ******************************************************************

      SUBROUTINE RUNOFF_HORTON()

      USE BIEF_DEF
      USE INTERFACE_TELEMAC2D, EX_RUNOFF_HORTON => RUNOFF_HORTON

      IMPLICIT NONE

!     Locals
      INTEGER(4) :: I
      REAL(8) :: WT_INFLTRTN

!     ******************************************************************
!     Users need to add declarations of followings, either as input/output argument or from exisiting subroutine that is used for declaration. For structures, users need to allocate memory before computation and deallocate memory after computation.
!     Users need to initialize all parameters and INFLTRTN_CMLT_CMPT.
!     ******************************************************************

!     NPOIN:                number of mesh nodes, (INTEGER*4)
!     INFLTRTN:             infiltration water, (TYPE(BIEF_OBJ))
!     INFLTRTN_RT_INTL:     parameter of initial infiltration rate, (TYPE(BIEF_OBJ))
!     H_PRE:                water depth of pre-computation, (TYPE(BIEF_OBJ))
!     HN:                   water depth of current time level, (TYPE(BIEF_OBJ))
!     DCY_COEFFCNT:         parameter of decay coefficient, (TYPE(BIEF_OBJ))
!     INFLTRTN_CMLT_CMPT:   cumulative excess infiltration water, (TYPE(BIEF_OBJ))
!     INFLTRTN_RT_FNL:      parameter of final infiltration rate, (TYPE(BIEF_OBJ))
!     DT:                   time step, (REAL*8)
!     RCVRY_COEFFCNT:       parameter of recovery coefficient, (TYPE(BIEF_OBJ))
!     OPTSOU:               type of sources, (INTEGER*4)
!     SMH:                  total mass source on right hand side of mass equation, (TYPE(BIEF_OBJ))
!     UNSV2D:               inverse of covering area of mesh nodes, (TYPE(BIEF_OBJ))

!     ******************************************************************

!     ******************************************************************
!     The excess infiltration is positive only when actual infiltration is greater than saturated hydraulic conduction. It is suggested that the part percolating into deeper soil layer is always a part of actual infiltration, and it is not related to recovery of infiltration capacity. When actual infiltration is less than or equal to saturated hydraulic conduction, the excess infiltration is set to zero, and cumulative excess infiltration does not vary.
!     As long as there is water, infiltration capacity does not recover, even if water amount is really small and actual infiltration is less than saturated hydraulic conduction.
!     If there is no water, there is no actual infiltration, and infiltration capacity recovers by decreasing cumulative excess infiltration.
!     ******************************************************************

      DO I = 1, NPOIN ! Node by node
         INFLTRTN%R(I) = 0.D0

         IF (INFLTRTN_RT_INTL%R(I) .LE. 0.D0) CYCLE

         IF ((H_PRE%R(I) .GT. 0.D0) .AND. (H_PRE%R(I) .GE. HN%R(I)))
     &   THEN
            WT_INFLTRTN = H_PRE%R(I)
         ELSEIF ((HN%R(I) .GT. 0.D0) .AND. (H_PRE%R(I) .LT. HN%R(I)))
     &   THEN
            IF (H_PRE%R(I) .LT. 0.D0) THEN
               WT_INFLTRTN =
     &         (0.5D0 * HN%R(I)) * (HN%R(I) / (HN%R(I) - H_PRE%R(I)))
            ELSE
               WT_INFLTRTN = 0.5D0 * (HN%R(I) + H_PRE%R(I))
            END IF
         ELSE
            WT_INFLTRTN = 0.D0
            GOTO 100
         END IF

         INFLTRTN%R(I) = MAX((INFLTRTN_RT_INTL%R(I)
     &   - (DCY_COEFFCNT%R(I) * INFLTRTN_CMLT_CMPT%R(I))),
     &   INFLTRTN_RT_FNL%R(I))
     &   * DT

         INFLTRTN%R(I) = MIN(INFLTRTN%R(I), WT_INFLTRTN)

100      IF ((WT_INFLTRTN .GT. 0.D0) .AND.
     &   (INFLTRTN%R(I) .GT. (INFLTRTN_RT_FNL%R(I) * DT))) THEN
            INFLTRTN_CMLT_CMPT%R(I) = INFLTRTN_CMLT_CMPT%R(I) +
     &      (INFLTRTN%R(I) - (INFLTRTN_RT_FNL%R(I) * DT))

         ELSEIF ((WT_INFLTRTN .LE. 0.D0) .AND.
     &   (INFLTRTN_CMLT_CMPT%R(I) .GT. 0.D0) .AND.
     &   (RCVRY_COEFFCNT%R(I) .GT. 0.D0)) THEN
            INFLTRTN_CMLT_CMPT%R(I) = INFLTRTN_CMLT_CMPT%R(I) *
     &      (EXP(- (RCVRY_COEFFCNT%R(I) * DT)))
         END IF

         IF (INFLTRTN%R(I) .GT. 0.D0) THEN
            IF (OPTSOU .EQ. 1) THEN
               SMH%R(I) = SMH%R(I) - (INFLTRTN%R(I) / DT)
            ELSE
               SMH%R(I) = SMH%R(I) - (INFLTRTN%R(I) / UNSV2D%R(I) / DT)
            END IF
         END IF
      END DO ! Node by node

      RETURN

      END SUBROUTINE