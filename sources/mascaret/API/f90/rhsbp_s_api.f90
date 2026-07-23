! obtain width, perimeter, flow area, and hydraulic radius
subroutine RHSBP_S_API(B1, B2, BS, P1, P2, S1, S2, R1, R2, &
                       Section, Z, ZREF, IDT, XDT, Profil, PROF, &
                       UniteListing, Erreur)

! declarations ---------------------------------------------------------------------------------------------------------------------
   use M_PRECISION
   use M_MESSAGE_C
   use M_PARAMETRE_C
   use M_PROFIL_T
   use M_PROFIL_PLAN_T
   use M_ERREUR_T
   use M_TRAITER_ERREUR_I

   implicit none

   ! outputs
   real(8),                     intent(out) :: B1, B2, BS, P1, P2, S1, S2, R1, R2
   type(ERREUR_T),            intent(inout) :: Erreur

   ! inputs
   type(PROFIL_T), dimension(:), intent(in) :: Profil
   type(PROFIL_PLAN_T),          intent(in) :: Prof
   integer(4),                   intent(in) :: Section
   real(8),                      intent(in) :: Z, ZREF
   real(8),        dimension(:), intent(in) :: XDT
   integer(4),     dimension(:), intent(in) :: IDT
   integer(4),                   intent(in) :: UniteListing

   ! locals
   integer(4)                               :: I, IP1, J, K
   real(8)                                  :: PAS, XD, Y, YD
   integer(4)                               :: NBPAS
   real(8),                    dimension(2) :: FB1, FB2, FBS, FP1, FP2, FS1, FS2

   intrinsic INT

! subroutine codes -----------------------------------------------------------------------------------------------------------------
   Erreur%Numero = 0

   S1  = W0
   S2  = W0
   B1  = W0
   B2  = W0
   BS  = W0
   R1  = W0
   R2  = W0
   FS1 = W0
   Y   = Z - ZREF

   if( Y <= EPS3 ) then
      Erreur%Numero = 1
      Erreur%ft   = err_1
      Erreur%ft_c = err_1c
      call TRAITER_ERREUR (Erreur, Section)
      return
   end if

   XD  = XDT(Section)
   I   = IDT(Section)
   IP1 = I + 1
   if( XD <= EPS6 ) then
      IP1 = I
   end if

   PAS   = Profil(I)%Pas + ( Profil(IP1)%Pas - Profil(I)%Pas ) * XD
   NBPAS = Profil(I)%NbPas

   K = INT( Y / PAS ) + 1
   if( K >= NBPAS ) then
      K = NBPAS - 1
      if (UniteListing>0) then
         write (UniteListing,10000) Section, Y, NBPAS, PAS
      endif
   end if

   YD = Y - ( K - 1 ) * PAS

   do J = 1 , 2
      FB1(J) = Prof%B1(I,K) + ( Prof%B1(IP1,K) - Prof%B1(I,K) ) * XD
      FB2(J) = Prof%B2(I,K) + ( Prof%B2(IP1,K) - Prof%B2(I,K) ) * XD
      FP1(J) = Prof%P1(I,K) + ( Prof%P1(IP1,K) - Prof%P1(I,K) ) * XD
      FP2(J) = Prof%P2(I,K) + ( Prof%P2(IP1,K) - Prof%P2(I,K) ) * XD
      FS1(J) = Prof%S1(I,K) + ( Prof%S1(IP1,K) - Prof%S1(I,K) ) * XD
      FS2(J) = Prof%S2(I,K) + ( Prof%S2(IP1,K) - Prof%S2(I,K) ) * XD
      FBS(J) = Prof%BS(I,K) + ( Prof%BS(IP1,K) - Prof%BS(I,K) ) * XD
      K = K + 1
   end do

   B1 = FB1(1) + ( FB1(2) - FB1(1) ) * YD / PAS
   B2 = FB2(1) + ( FB2(2) - FB2(1) ) * YD / PAS
   P1 = FP1(1) + ( FP1(2) - FP1(1) ) * YD / PAS
   P2 = FP2(1) + ( FP2(2) - FP2(1) ) * YD / PAS
   S1 = FS1(1) + ( FS1(2) - FS1(1) ) * YD / PAS
   S2 = FS2(1) + ( FS2(2) - FS2(1) ) * YD / PAS
   BS = FBS(1) + ( FBS(2) - FBS(1) ) * YD / PAS
   R1 = S1 / P1

   if( BS < EPS3 ) then
      BS = 0._DOUBLE
   end if

   if( P2 > EPS3 ) then
      R2 = S2 / P2
   else
      R2 = 0._DOUBLE
      S2 = 0._DOUBLE
      P2 = 0._DOUBLE
      B2 = 0._DOUBLE
   end if

   return

   10000 format (                                                               &
        '<< ATTENTION >>',                                                      &
        'Dans la section de calcul n0 ',i4,/,                                   &
        'Tirant d''eau = ',g15.7,' depassant la hauteur du profil',/,           &
        'Augmenter le nombre de pas de planimetrage : ',i4,' ou le pas : ',f7.2)

end subroutine RHSBP_S_API