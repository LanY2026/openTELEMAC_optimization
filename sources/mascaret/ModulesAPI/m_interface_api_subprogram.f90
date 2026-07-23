! interfaces of api subprograms called by upper-level api programs
MODULE M_INTERFACE_API_SUBPROGRAM

   INTERFACE
   SUBROUTINE RHSBP_S_API(B1, B2, BS, P1, P2, S1, S2, R1, R2, &
                          Section, Z, ZREF, IDT, XDT, Profil, PROF, &
                          UniteListing, Erreur)
   use M_ERREUR_T
   use M_PROFIL_T
   use M_PROFIL_PLAN_T
   implicit none
   real(8),                     intent(out) :: B1, B2, BS, P1, P2, S1, S2, R1, R2
   type(ERREUR_T),            intent(inout) :: Erreur
   type(PROFIL_T), dimension(:), intent(in) :: Profil
   type(PROFIL_PLAN_T),          intent(in) :: Prof
   integer(4),                   intent(in) :: Section
   real(8),                      intent(in) :: Z, ZREF
   real(8),        dimension(:), intent(in) :: XDT
   integer(4),     dimension(:), intent(in) :: IDT
   integer(4),                   intent(in) :: UniteListing
   END SUBROUTINE RHSBP_S_API
   END INTERFACE

   INTERFACE
   SUBROUTINE REPAR_API(DEB, BETA, &
                        S1, S2, RH1, RH2, &
                        P1, P2, CF1, CF2, ModeleLit, LoiFrottement, NomProfil, &
                        Erreur)
   use M_ERREUR_T
   implicit none
   real(8),        intent(out)   :: DEB, BETA
   real(8),        intent(inout) :: S1, S2, RH1, RH2
   type(ERREUR_T), intent(inout) :: Erreur
   real(8),        intent(in)    :: P1, P2, CF1, CF2
   integer(4),     intent(in)    :: ModeleLit, LoiFrottement
   Character(30),  intent(in)    :: NomProfil
   END SUBROUTINE REPAR_API
   END INTERFACE

END MODULE M_INTERFACE_API_SUBPROGRAM