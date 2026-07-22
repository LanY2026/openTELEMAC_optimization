!     Authored by Ziyi Huang, December 2025
!     Computation of infiltration water node by node, in coupling modeling of TELEMAC-2D (version tag v9p0r0) and SCS curve number infiltration model. SMH is revised if there is infiltration water.

!     ******************************************************************
!     This subroutine should be called by a subroutine that sets SMH, which generally is PROSOU.
!     Execution of pre-computation is required before this subroutine is called.
!     ******************************************************************

      SUBROUTINE RUNOFF_SCS_CN()

      USE BIEF_DEF
      USE INTERFACE_TELEMAC2D, EX_RUNOFF_SCS_CN => RUNOFF_SCS_CN

      IMPLICIT NONE

!     locals
      INTEGER(4) :: I
      REAL(8) :: INFLTRTN_CMLT_NEW

!     ******************************************************************
!     Users need to add declarations of followings, either as input/output argument or from exisiting subroutine that is used for declaration. For structures, users need to allocate memory before computation and deallocate memory after computation.
!     Users need to initialize all parameters and WT_INFLTRTN_CMLT, INFLTRTN_CMLT_CMPT.
!     ******************************************************************

!     NPOIN:                number of mesh nodes, (INTEGER*4)
!     INFLTRTN:             infiltration water, (TYPE(BIEF_OBJ))
!     INFLTRTN_CMLT_MAX:    parameter of maximum of cumulative infiltration water of increasing part, (TYPE(BIEF_OBJ))
!     H_PRE:                water depth of pre-computation, (TYPE(BIEF_OBJ))
!     HN:                   water depth of current time level, (TYPE(BIEF_OBJ))
!     WT_INFLTRTN_CMLT:     cumulation of increasing part of water source for infiltration, (TYPE(BIEF_OBJ))
!     INFLTRTN_CMLT_CMPT:   cumulative infiltration water of increasing part, (TYPE(BIEF_OBJ))
!     KS:                   parameter of saturated hydraulic conductivity, (TYPE(BIEF_OBJ))
!     DT:                   time step, (REAL*8)
!     OPTSOU:               type of sources, (INTEGER*4)
!     SMH:                  total mass source on right hand side of mass equation, (TYPE(BIEF_OBJ))
!     UNSV2D:               inverse of covering area of mesh nodes, (TYPE(BIEF_OBJ))

!     ******************************************************************

!     ******************************************************************
!     For increasing part of water, use SCS curve number model. For existing part of water, infiltration is limited by saturated hydraulic conductivity.
!     ******************************************************************

      DO I = 1, NPOIN ! Node by node
         INFLTRTN%R(I) = 0.D0

         IF ((INFLTRTN_CMLT_MAX%R(I) .GT. 0.D0) .AND.
     &   (H_PRE%R(I) .GT. 0.D0) .AND.
     &   (H_PRE%R(I) .GT. HN%R(I))) THEN
            WT_INFLTRTN_CMLT%R(I) =
     &      WT_INFLTRTN_CMLT%R(I) + (H_PRE%R(I) - MAX(HN%R(I), 0.D0))
            INFLTRTN_CMLT_NEW =
     &      (WT_INFLTRTN_CMLT%R(I) * INFLTRTN_CMLT_MAX%R(I)) /
     &      (WT_INFLTRTN_CMLT%R(I) + INFLTRTN_CMLT_MAX%R(I))

            INFLTRTN%R(I) = INFLTRTN%R(I) +
     &      (INFLTRTN_CMLT_NEW - INFLTRTN_CMLT_CMPT%R(I))

            INFLTRTN_CMLT_CMPT%R(I) = INFLTRTN_CMLT_NEW
         END IF

         IF ((KS%R(I) .GT. 0.D0) .AND. (HN%R(I) .GT. 0.D0)) THEN
            IF (H_PRE%R(I) .LT. HN%R(I)) THEN
               IF (H_PRE%R(I) .LT. 0.D0) THEN
                  INFLTRTN%R(I) = INFLTRTN%R(I) +
     &            MIN(((0.5D0 * HN%R(I)) *
     &            (HN%R(I) / (HN%R(I) - H_PRE%R(I)))),
     &            (KS%R(I) * DT))

               ELSE
                  INFLTRTN%R(I) = INFLTRTN%R(I) +
     &            MIN((0.5D0 * (HN%R(I) + H_PRE%R(I))), (KS%R(I) * DT))
               END IF

            ELSE
               INFLTRTN%R(I) =
     &         INFLTRTN%R(I) + MIN(HN%R(I), (KS%R(I) * DT))
            END IF
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