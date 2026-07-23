! obtain conveyance and beta, and revise flow-area and hydraulic radius
subroutine REPAR_API(DEB, BETA, &
                     S1, S2, RH1, RH2, &
                     P1, P2, CF1, CF2, ModeleLit, LoiFrottement, NomProfil, &
                     Erreur)

! declarations ---------------------------------------------------------------------------------------------------------------------
   use M_PRECISION
   use M_PARAMETRE_C
   use M_MESSAGE_C
   use M_CONSTANTES_CALCUL_C
   use M_ERREUR_T
   use M_DEBITANCE_S
   use M_TRAITER_ERREUR_I

   implicit none

   ! outputs
   real(8),        intent(out)   :: DEB, BETA
   real(8),        intent(inout) :: S1, S2, RH1, RH2
   type(ERREUR_T), intent(inout) :: Erreur

   ! inputs
   real(8),        intent(in)    :: P1, P2, CF1, CF2
   integer(4),     intent(in)    :: ModeleLit, LoiFrottement
   Character(30),  intent(in)    :: NomProfil

   ! locals
   real(8)                       :: st1_temp
   real(8)                       :: RH
   real(8)                       :: A, A0, DEB1, DEB2, ETA, FP1, FP2, FS1, FS2, R0, S, STEQUI, USETA, VALOP
   real(8),            PARAMETER :: PUT = 0.3_DOUBLE

   intrinsic DCOS, DSQRT

! subroutine codes -----------------------------------------------------------------------------------------------------------------
   Erreur%Numero = 0

   if( P1 <= EPS6 ) then
      Erreur%Numero = 40
      Erreur%ft     = err_40
      Erreur%ft_c   = err_40c
      call TRAITER_ERREUR_1C1R( Erreur , NomProfil, P1 )
      return
   endif

   if( S1 <= EPS6 ) then
      Erreur%Numero = 41
      Erreur%ft     = err_41
      Erreur%ft_c   = err_41c
      call TRAITER_ERREUR( Erreur , S1 )
      return
   endif

   RH   = ( S1 + S2 ) / ( P1 + P2 )

   label_FOND_BERGE : if( ModeleLit == MODELE_LIT_FOND_BERGE ) then

      S      = S1 + S2
      FP1    = P1 / ( CF1**W32 )
      FP2    = P2 / ( CF2**W32 )
      FS1    = S * ( FP1 / ( FP1 + FP2 ) )
      FS2    = S * ( FP2 / ( FP1 + FP2 ) )
      USETA  = FS2 / FS1
      STEQUI = ( ( P1 + P2 ) / ( FP1 + FP2 ) )**W23

      DEB  = STEQUI * S * RH**W23
      BETA = 1._DOUBLE

      S1 = FS1
      S2 = FS2
      RH1 = S1 / P1
      if( abs(P2).GT.EPS6 ) then
         RH2 = S2 / P2
      else
         RH2 = 0._DOUBLE
      end if

   else label_FOND_BERGE

      A = 1._DOUBLE

      if( ModeleLit == MODELE_LIT_DEBORD ) then
         A0 = W09 * ( CF2 / CF1 )**W16
         R0 = RH2 / RH1
         if( R0 >= PUT ) then
            A = A0
         else
            A = ( ( 1._DOUBLE - A0 ) * DCOS( PI * R0 / PUT ) + 1._DOUBLE + A0 ) / 2._DOUBLE
         end if
      end if

      call DEBITANCE_S(DEB1, st1_temp, &
                       RH1, S1, LoiFrottement, CF1, &
                       Erreur)
      if( Erreur%Numero /= 0 ) then
         return
      endif

      DEB1 = DEB1 * A
      VALOP = S2**2 + S1 * S2 * ( 1._DOUBLE - A * A )
      IF(VALOP.LT.0.D0) THEN
         Erreur%Numero = 705
         Erreur%ft     = err_705
         Erreur%ft_c   = err_705c
         call TRAITER_ERREUR( Erreur , VALOP )
         return
      ENDIF
      DEB2 = CF2 * DSQRT( VALOP ) * RH2**W23
      DEB  = DEB1 + DEB2

      if( S2 <= ( S1 * EPS4 ) ) then
         BETA = 1._DOUBLE
      else
         ETA  = DEB1 / DEB2
         BETA = ( ETA**2 / S1 + 1._DOUBLE / S2 ) * ( S1 + S2 ) / ( 1._DOUBLE + ETA )**2
      end if

   end if label_FOND_BERGE

end subroutine REPAR_API