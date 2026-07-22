!     Authored by Ziyi Huang, February 2026
!     Computation of infiltration water node by node, in coupling modeling of TELEMAC-2D (version tag v9p0r0) and Green-Ampt infiltration model. SMH is revised if there is infiltration water.

!     ******************************************************************
!     This subroutine should be called by a subroutine that sets SMH, which generally is PROSOU.
!     Execution of pre-computation is required before this subroutine is called.
!     ******************************************************************

      SUBROUTINE RUNOFF_GREENAMPT()

      USE BIEF
      USE BIEF_DEF
      USE INTERFACE_TELEMAC2D, EX_RUNOFF_GREENAMPT => RUNOFF_GREENAMPT

      IMPLICIT NONE

!     Locals
      INTEGER(4) :: I
      REAL(8) :: WT_INFLTRTN, INFLTRTN_CMLT_STRTN

!     ******************************************************************
!     Users need to add declarations of followings, either as input/output argument or from exisiting subroutine that is used for declaration. For structures, users need to allocate memory before computation and deallocate memory after computation.
!     Users need to initialize all parameters and STRTN, TM_RCVRY_RMN, MSTR_DFCT, MSTR_DFCT_INTL, INFLTRTN_CMLT_CMPT.
!     ******************************************************************

!     NPOIN:                number of mesh nodes, (INTEGER*4)
!     INFLTRTN:             infiltration water, (TYPE(BIEF_OBJ))
!     KS:                   parameter of saturated hydraulic conductivity, (TYPE(BIEF_OBJ))
!     H_PRE:                water depth of pre-computation, (TYPE(BIEF_OBJ))
!     HN:                   water depth of current time level, (TYPE(BIEF_OBJ))
!     STRTN:                there is saturated layer or not, (TYPE(BIEF_OBJ))
!     TM_RCVRY_RMN:         remaining time of recovery, (TYPE(BIEF_OBJ))
!     DT:                   time step, (REAL*8)
!     MSTR_DFCT:            moisture content deficit to saturated state, (TYPE(BIEF_OBJ))
!     MSTR_DFCT_MAX:        parameter of maximum of subsurface moisture content deficit to saturated state, (TYPE(BIEF_OBJ))
!     RCVRY_COEFFCNT:       parameter of recovery coefficient, (TYPE(BIEF_OBJ))
!     MSTR_DFCT_INTL:       initial subsurface moisture content deficit to saturated state, (TYPE(BIEF_OBJ))
!     INFLTRTN_CMLT_CMPT:   cumulative infiltration water, (TYPE(BIEF_OBJ))
!     DPTH_UPPR_ZN:         parameter of thickness of upper layer of subsurface that influences infiltration capacity, (TYPE(BIEF_OBJ))
!     TM_RCVRY:             parameter of recovery period for a new wet process, (TYPE(BIEF_OBJ))
!     SCTN_HD:              parameter of capillary suction head between saturated layer and unsaturated layer, (TYPE(BIEF_OBJ))
!     OPTSOU:               type of sources, (INTEGER*4)
!     SMH:                  total mass source on right hand side of mass equation, (TYPE(BIEF_OBJ))
!     UNSV2D:               inverse of covering area of mesh nodes, (TYPE(BIEF_OBJ))

!     ******************************************************************

!     ******************************************************************
!     Remaining time of recovery of infiltration capacity being zero means that infiltration capacity is fully recovered. Remaining time of recovery keeps decreasing (infiltration capacity keeps recovering) as long as there is no saturation. When time of recovery decreases to zero (infiltration capacity fully recovered) and there is no surface water, a new wet process starts for which cumulative infiltration is reset to zero and initial moisture deficit is reset. Whenever there is saturated layer or surface water exceeds saturated hydraulic conduction, remaining time of recovery increases back to recovery period for a new wet process.
!     The variable moisture content deficit to saturated state should not be less than zero and should not be greater than the constant maximum.
!     If saturating occurs in a time-step, infiltration solved in terms of Green-Ampt method does not exceed water source for infiltration in this time-step theorecically. Otherwise it must be numerical error.
!     ******************************************************************

      DO I = 1, NPOIN ! Node by node
         INFLTRTN%R(I) = 0.D0

         IF (KS%R(I) .LE. 0.D0) CYCLE

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
         END IF

         IF (STRTN%I(I) .EQ. 0) THEN ! There is saturated layer or not
            IF (WT_INFLTRTN .LE. 0.D0) THEN ! Water amount
               IF (TM_RCVRY_RMN%R(I) .GT. 0.D0)
     &         TM_RCVRY_RMN%R(I) = TM_RCVRY_RMN%R(I) - DT

               IF ((MSTR_DFCT%R(I) .LT. MSTR_DFCT_MAX%R(I))
     &         .AND. (RCVRY_COEFFCNT%R(I) .GT. 0.D0))
     &         MSTR_DFCT%R(I) =
     &         MIN((MSTR_DFCT%R(I) +
     &         (RCVRY_COEFFCNT%R(I) * DT * MSTR_DFCT_MAX%R(I))),
     &         MSTR_DFCT_MAX%R(I))

               IF (TM_RCVRY_RMN%R(I) .LE. 0.D0) THEN
                  IF (MSTR_DFCT_INTL%R(I) .NE. MSTR_DFCT%R(I))
     &            MSTR_DFCT_INTL%R(I) = MSTR_DFCT%R(I)

                  IF (INFLTRTN_CMLT_CMPT%R(I) .GT. 0.D0)
     &            INFLTRTN_CMLT_CMPT%R(I) = 0.D0

               ELSE
                  IF ((INFLTRTN_CMLT_CMPT%R(I) .GT. 0.D0) .AND.
     &            (RCVRY_COEFFCNT%R(I) .GT. 0.D0))
     &            INFLTRTN_CMLT_CMPT%R(I) =
     &            MAX((INFLTRTN_CMLT_CMPT%R(I) -
     &            (RCVRY_COEFFCNT%R(I) * DT *
     &            MSTR_DFCT_MAX%R(I) * DPTH_UPPR_ZN%R(I))),
     &            0.D0)
               END IF

            ELSEIF ((WT_INFLTRTN .GT. 0.D0) .AND.
     &      (WT_INFLTRTN .LE. (KS%R(I) * DT))) THEN ! Water amount
               IF (TM_RCVRY_RMN%R(I) .GT. 0.D0)
     &         TM_RCVRY_RMN%R(I) = TM_RCVRY_RMN%R(I) - DT

               INFLTRTN%R(I) = WT_INFLTRTN

               INFLTRTN_CMLT_CMPT%R(I) =
     &         INFLTRTN_CMLT_CMPT%R(I) + INFLTRTN%R(I)

               IF (MSTR_DFCT%R(I) .GT. 0.D0)
     &         MSTR_DFCT%R(I) =
     &         MAX((MSTR_DFCT%R(I) -
     &         (INFLTRTN%R(I) / DPTH_UPPR_ZN%R(I))),
     &         0.D0)

            ELSEIF (WT_INFLTRTN .GT. (KS%R(I) * DT)) THEN ! Water amount
               IF (TM_RCVRY_RMN%R(I) .NE. TM_RCVRY%R(I))
     &         TM_RCVRY_RMN%R(I) = TM_RCVRY%R(I)

               INFLTRTN_CMLT_STRTN =
     &         (KS%R(I) * MSTR_DFCT_INTL%R(I) * SCTN_HD%R(I)) /
     &         ((WT_INFLTRTN / DT) - KS%R(I))

               IF (INFLTRTN_CMLT_STRTN .GT.
     &         (INFLTRTN_CMLT_CMPT%R(I) + WT_INFLTRTN)) THEN ! Infiltration for saturation
                  INFLTRTN%R(I) = WT_INFLTRTN

                  INFLTRTN_CMLT_CMPT%R(I) =
     &            INFLTRTN_CMLT_CMPT%R(I) + INFLTRTN%R(I)

               ELSEIF ((INFLTRTN_CMLT_STRTN .LE.
     &         (INFLTRTN_CMLT_CMPT%R(I) + WT_INFLTRTN)) .AND.
     &         (INFLTRTN_CMLT_STRTN .GE.
     &         INFLTRTN_CMLT_CMPT%R(I))) THEN ! Infiltration for saturation
                  STRTN%I(I) = 1

                  INFLTRTN%R(I) = WT_INFLTRTN

                  INFLTRTN_CMLT_CMPT%R(I) =
     &            INFLTRTN_CMLT_CMPT%R(I) + INFLTRTN%R(I)

               ELSEIF (INFLTRTN_CMLT_STRTN .LT.
     &         INFLTRTN_CMLT_CMPT%R(I)) THEN ! Infiltration for saturation
                  STRTN%I(I) = 1

                  INFLTRTN%R(I) =
     &            (KS%R(I) *
     &            (1.D0 +
     &            ((MSTR_DFCT_INTL%R(I) * SCTN_HD%R(I)) /
     &            INFLTRTN_CMLT_CMPT%R(I)))) *
     &            DT

                  IF (INFLTRTN%R(I) .GT. WT_INFLTRTN) THEN
                     WRITE(LU, 100) I
                     CALL PLANTE(1)
                  END IF

                  INFLTRTN_CMLT_CMPT%R(I) =
     &            INFLTRTN_CMLT_CMPT%R(I) + INFLTRTN%R(I)
               END IF ! Infiltration for saturation

               IF (MSTR_DFCT%R(I) .GT. 0.D0)
     &         MSTR_DFCT%R(I) =
     &         MAX((MSTR_DFCT%R(I) -
     &         (INFLTRTN%R(I) / DPTH_UPPR_ZN%R(I))),
     &         0.D0)
            END IF ! Water amount

         ELSEIF (STRTN%I(I) .EQ. 1) THEN ! There is saturated layer or not
            IF (TM_RCVRY_RMN%R(I) .NE. TM_RCVRY%R(I))
     &      TM_RCVRY_RMN%R(I) = TM_RCVRY%R(I)

            IF (WT_INFLTRTN .LE. 0.D0) THEN ! Water amount
               STRTN%I(I) = 0

            ELSE ! Water amount
               INFLTRTN%R(I) =
     &         (KS%R(I) *
     &         (1.D0 +
     &         ((MSTR_DFCT_INTL%R(I) * MAX(HN%R(I), 0.D0)) /
     &         INFLTRTN_CMLT_CMPT%R(I)) +
     &         ((MSTR_DFCT_INTL%R(I) * SCTN_HD%R(I)) /
     &         INFLTRTN_CMLT_CMPT%R(I)))) *
     &         DT

               IF (INFLTRTN%R(I) .GT. WT_INFLTRTN) THEN
                  STRTN%I(I) = 0

                  INFLTRTN%R(I) = WT_INFLTRTN
               END IF

               INFLTRTN_CMLT_CMPT%R(I) =
     &         INFLTRTN_CMLT_CMPT%R(I) + INFLTRTN%R(I)

               IF (MSTR_DFCT%R(I) .GT. 0.D0)
     &         MSTR_DFCT%R(I) =
     &         MAX((MSTR_DFCT%R(I) -
     &         (INFLTRTN%R(I) / DPTH_UPPR_ZN%R(I))),
     &         0.D0)
            END IF ! Water amount
         END IF ! There is saturated layer or not

         IF (INFLTRTN%R(I) .GT. 0.D0) THEN
            IF (OPTSOU .EQ. 1) THEN
               SMH%R(I) = SMH%R(I) - (INFLTRTN%R(I) / DT)
            ELSE
               SMH%R(I) = SMH%R(I) - (INFLTRTN%R(I) / UNSV2D%R(I) / DT)
            END IF
         END IF
      END DO ! Node by node

100   FORMAT('TELEMAC2D - RUNOFF_GREENAMPT: ',
     &'numerical solution of infiltration by Green-Ampt method ',
     &'is too large at node No. ' I10, '.')

      RETURN

      END SUBROUTINE