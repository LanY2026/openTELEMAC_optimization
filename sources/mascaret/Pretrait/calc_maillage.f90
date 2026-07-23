!== Copyright (C) 2000-2022 EDF-CEREMA ==
!
!   This file is part of MASCARET.
!
!   MASCARET is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   MASCARET is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with MASCARET.  If not, see <http://www.gnu.org/licenses/>
!

subroutine CALC_MAILLAGE ( &
     X                   , & ! Tableau des abscisses
     TypeMaillage        , & ! Type de calcul du maillage
     FichierMaillage     , & ! Fichier du maillage
     FichierSauveMaillage, & ! Fichier de sauvegarde du maillage
     Profil              , & ! Profils geometriques
     ProfDebBief         , & ! Premier profil d'un bief
     ProfFinBief         , & ! Dernier profil d'un bief
     AbscRelExtDebBief   , & ! Abscisse rel de l'extremite debut du bief
     AbscRelExtFinBief   , & ! Abscisse rel de l'extremite debut du bief
     impression_geo      , & ! Flag d'impression de la geometrie
     UniteListing        , & ! Unite logique fichier listing
     unitNum             , & ! Unite logique .xcas
     Erreur                & ! Erreur
                         )
! *********************************************************************
! PROGICIEL : MASCARET         A. LEBOSSE
!                              S. MANDELKERN
!                              F. ZAOUI
!
! VERSION : V8P4R0                EDF-CEREMA
! *********************************************************************

   !========================= Declarations ===========================
   use M_PRECISION
   use M_PARAMETRE_C
   use M_ERREUR_T            ! Type ERREUR_T
   use M_FICHIER_T            ! Type FICHIER_T
   use M_MAILLE_T            ! Types MAILLE_E_T et MAILLE_R_T
   use M_PROFIL_T            ! Type PROFIL_T
   use M_MESSAGE_C           ! Messages d'erreur
   use M_CONSTANTES_CALCUL_C ! Constantes num, phys et info
   use M_LEC_MAILLAGE_I      ! Lecture du maillage
   use M_MAILLER_I           ! Interface de sous-programme
   use M_TRAITER_ERREUR_I    ! Traitement de l'errreur
   use M_ABS_ABS_S           ! Calcul de l'abscisse absolue
   use M_XCAS_S

   implicit none

   type SECTION_REL_T
      sequence
      integer      :: Branche     ! Numero de branche
      real(DOUBLE) :: AbscisseRel ! Abscisse relative
   end type SECTION_REL_T

   ! Arguments
   real(DOUBLE)      , dimension(:), pointer       :: X
   integer                         , intent(  out) :: TypeMaillage
   type(FICHIER_T)                 , intent(inout) :: FichierMaillage
   type(FICHIER_T)                 , intent(inout) :: FichierSauveMaillage
   type(PROFIL_T)    , dimension(:), intent(in   ) :: Profil
   integer           , dimension(:), intent(in   ) :: ProfDebBief
   integer           , dimension(:), intent(in   ) :: ProfFinBief
   logical                         , intent(in   ) :: impression_geo
   integer                         , intent(in   ) :: UniteListing
   real(DOUBLE)      , dimension(:), intent(in   ) :: AbscRelExtDebBief
   real(DOUBLE)      , dimension(:), intent(in   ) :: AbscRelExtFinBief
   integer, intent(in)                             :: unitNum
   type(ERREUR_T)                  , intent(inout) :: Erreur
   ! Variables locales
   type(SECTION_REL_T) :: section
   integer             :: mode_saisie_maillage
   integer             :: nb_section
   integer             :: nb_maille
   logical             :: sauvegarde_maillage
   type(MAILLE_R_T), dimension(:), pointer :: maille_r => null()
   type(MAILLE_E_T), dimension(:), pointer :: maille_e => null()
   integer :: retour              ! code de retour des fonctions
                                  ! intrinseques
   integer :: k                   ! compteur sur les mailles
   integer :: imail               ! compteur sur les mailles
   integer, allocatable :: itab(:)
   real(double), allocatable :: rtab(:)
   character(len=256)  :: pathNode
   character(len=8192) :: line
   !character(132) :: !arbredappel_old

   !========================= Instructions ===========================
   ! INITIALISATION
   ! --------------
   Erreur%Numero = 0
   retour = 0
   k      = 0
   !arbredappel_old = trim(!Erreur%arbredappel)
   !Erreur%arbredappel = trim(!Erreur%arbredappel)//'=>CALC_MAILLAGE'

   if (UniteListing >0) write(UniteListing,10465)

   ! Mode de saisie du maillage
   !---------------------------
   pathNode = 'parametresPlanimetrageMaillage/maillage/modeSaisie'
   line = xcasReader(unitNum, pathNode)
   read(unit=line, fmt=*) mode_saisie_maillage

   if( mode_saisie_maillage /= SAISIE_PAR_FICHIER .and. &
      mode_saisie_maillage /= SAISIE_PAR_CLAVIER) then
      Erreur%Numero = 305
      Erreur%ft     = err_305
      Erreur%ft_c   = err_305c
      call TRAITER_ERREUR( Erreur , 'Mode de saisie du maillage' , '1 ou 2' )
      return
   end if

   if( mode_saisie_maillage == SAISIE_PAR_CLAVIER ) then
      if (UniteListing >0) write(UniteListing,10480) 'PAR CLAVIER'
      ! Methode de calcul du maillage
      pathNode = 'parametresPlanimetrageMaillage/methodeMaillage'
      line = xcasReader(unitNum, pathNode)
      if(len(trim(line)).eq.0) then
        print*,"Parse error => methodeMaillage"
        call xerror(Erreur)
        return
      endif
      read(unit=line, fmt=*) TypeMaillage
      if( TypeMaillage < 1 .or. TypeMaillage > TYPE_MAILLAGE_NB_MAX ) then
         Erreur%Numero = 371
         Erreur%ft   = err_371
         Erreur%ft_c = err_371c
         call TRAITER_ERREUR( Erreur , TYPE_MAILLAGE_NB_MAX )
         return
      end if

      select case(TypeMaillage)
         case( TYPE_MAILLAGE_PROFIL )
            if (UniteListing >0) write(UniteListing,10470) 'AUX PROFILS'
         case( TYPE_MAILLAGE_SERIE )
            if (UniteListing >0) write(UniteListing,10470) 'PAR SERIES'
         case( TYPE_MAILLAGE_INDIVIDUELLE )
            if (UniteListing >0) write(UniteListing,10470) 'SECTION PAR SECTION'
         case( TYPE_MAILLAGE_PRECEDENT )
            if (UniteListing >0) write(UniteListing,10470) 'REPRISE DU MAILLAGE PRECEDENT'
         case(TYPE_MAILLAGE_SERIE_PROFIL)
            if (UniteListing >0) write(UniteListing,10470) 'PAR SERIE AUX PROFILS'
      end select

      select case( TypeMaillage )
         case( TYPE_MAILLAGE_SERIE )
            pathNode = 'parametresPlanimetrageMaillage/maillage/maillageClavier/nbZones'
            line = xcasReader(unitNum, pathNode)
            read(unit=line, fmt=*) nb_maille

            if( nb_maille <= 0 ) then
               Erreur%Numero = 306
               Erreur%ft     = err_306
               Erreur%ft_c   = err_306c
               call TRAITER_ERREUR( Erreur , 'Nombre de zones pour le calcul du maillage' )
               return
            end if

            if(.not.associated(maille_r)) allocate( maille_r(nb_maille) , STAT = retour )
            if( retour /= 0 ) then
               Erreur%Numero = 5
               Erreur%ft     = err_5
               Erreur%ft_c   = err_5c
               call TRAITER_ERREUR( Erreur , 'maille_r' )
               return
            end if

            if(.not.associated(maille_e)) allocate( maille_e(0) , STAT = retour )
            if( retour /= 0 ) then
               Erreur%Numero = 5
               Erreur%ft     = err_5
               Erreur%ft_c   = err_5c
               call TRAITER_ERREUR( Erreur , 'maille_e' )
               return
            end if

            pathNode = 'parametresPlanimetrageMaillage/maillage/maillageClavier/numBrancheZone'
            line = xcasReader(unitNum, pathNode)
            read(unit=line, fmt=*) maille_r%Branche

            pathNode = 'parametresPlanimetrageMaillage/maillage/maillageClavier/absDebutZone'
            line = xcasReader(unitNum, pathNode)
            read(unit=line, fmt=*) maille_r%AbscisseDeb

            pathNode = 'parametresPlanimetrageMaillage/maillage/maillageClavier/absFinZone'
            line = xcasReader(unitNum, pathNode)
            read(unit=line, fmt=*) maille_r%AbscisseFin

            pathNode = 'parametresPlanimetrageMaillage/maillage/maillageClavier/nbSectionZone'
            line = xcasReader(unitNum, pathNode)
            read(unit=line, fmt=*) maille_r%NbSection

            RCHID_1: DO K = 1, SIZE(PROFDEBBIEF)
               DO IMAIL = 1, NB_MAILLE
                  IF (MAILLE_R(IMAIL)%BRANCHE .EQ. K) GOTO 100
               ENDDO
               ERREUR%NUMERO = 800
               ERREUR%FT = ERR_800
               ERREUR%FT_C = ERR_800C
               CALL TRAITER_ERREUR(ERREUR, 'NUMERO DE BRANCHE DE ZONE', K)
               RETURN
100            DO IMAIL = 1, NB_MAILLE
                  IF ((MAILLE_R(IMAIL)%BRANCHE .EQ. K) .AND. (MAILLE_R(IMAIL)%ABSCISSEDEB .EQ. ABSCRELEXTDEBBIEF(K))) GOTO 200
               ENDDO
               ERREUR%NUMERO = 801
               ERREUR%FT = ERR_801
               ERREUR%FT_C = ERR_801C
               CALL TRAITER_ERREUR(ERREUR, 'ABSCISSE DE DEBUT DE ZONE', K)
               RETURN
200            DO IMAIL = 1, NB_MAILLE
                  IF ((MAILLE_R(IMAIL)%BRANCHE .EQ. K) .AND. (MAILLE_R(IMAIL)%ABSCISSEFIN .EQ. ABSCRELEXTFINBIEF(K))) CYCLE RCHID_1
               ENDDO
               ERREUR%NUMERO = 801
               ERREUR%FT = ERR_801
               ERREUR%FT_C = ERR_801C
               CALL TRAITER_ERREUR(ERREUR, 'ABSCISSE DE FIN DE ZONE', K)
               RETURN
            ENDDO RCHID_1

            do k = 1 , nb_maille

               IF ((MAILLE_R(K)%BRANCHE .LT. 1) .OR. (MAILLE_R(K)%BRANCHE .GT. SIZE(PROFDEBBIEF))) THEN
                  ERREUR%NUMERO = 802
                  ERREUR%FT = ERR_802
                  ERREUR%FT_C = ERR_802C
                  CALL TRAITER_ERREUR(ERREUR, 'NUMERO DE BRANCHE DE ZONE', K)
                  RETURN
               ENDIF

               if( maille_r(k)%AbscisseDeb < AbscRelExtDebBief(maille_r(k)%Branche ) .or. &
                   maille_r(k)%AbscisseDeb > AbscRelExtFinBief(maille_r(k)%Branche)) then
                  Erreur%Numero   = 334
                  Erreur%ft       = err_334
                  Erreur%ft_c     = err_334c
                  call TRAITER_ERREUR( Erreur , 'debut' , k , maille_r(k)%Branche )
                  return
               end if

               if( maille_r(k)%AbscisseFin < AbscRelExtDebBief(maille_r(k)%Branche) .or. &
                   maille_r(k)%AbscisseFin > AbscRelExtFinBief(maille_r(k)%Branche)) then
                  Erreur%Numero   = 334
                  Erreur%ft       = err_334
                  Erreur%ft_c     = err_334c
                  call TRAITER_ERREUR( Erreur , 'fin' , k , maille_r(k)%Branche )
                  return
               end if

               if( maille_r(k)%AbscisseFin <= maille_r(k)%AbscisseDeb ) then
                  Erreur%Numero = 314
                  Erreur%ft     = err_314
                  Erreur%ft_c   = err_314c
                  call TRAITER_ERREUR( Erreur , 'des zones de maillage' , k )
                  return
               end if

               if( maille_r(k)%NbSection < 0 ) then
                  Erreur%Numero = 374
                  Erreur%ft     = err_374
                  Erreur%ft_c   = err_374c
                  call TRAITER_ERREUR( Erreur , 'du nombre de sections d''une maille du maillage' , k )
                  return
               end if

               ! Passage en abscisses absolues
               !------------------------------
               maille_r(k)%AbscisseDeb = ABS_ABS_S ( &
                           maille_r(k)%Branche     , &
                           maille_r(k)%AbscisseDeb , &
                           Profil                  , &
                           ProfDebBief             , &
                           ProfFinBief             , &
                           Erreur                    &
                                                 )
               if( Erreur%Numero /= 0 ) then
                  return
               end if

               maille_r(k)%AbscisseFin = ABS_ABS_S  ( &
                            maille_r(k)%Branche     , &
                            maille_r(k)%AbscisseFin , &
                            Profil                  , &
                            ProfDebBief             , &
                            ProfFinBief             , &
                            Erreur                    &
                                                 )
               if( Erreur%Numero /= 0 ) then
                  return
               end if
            end do

            ! Controle de non chevauchement des mailles
            !------------------------------------------
            do imail = 2 , nb_maille
               IF (MAILLE_R(IMAIL)%BRANCHE .LT. MAILLE_R(IMAIL - 1)%BRANCHE) THEN
                  ERREUR%NUMERO = 803
                  ERREUR%FT = ERR_803
                  ERREUR%FT_C = ERR_803C
                  CALL TRAITER_ERREUR(ERREUR, 'NUMERO DE BRANCHE DE ZONE', (IMAIL - 1), IMAIL)
                  RETURN
               ENDIF
               if( maille_r(imail)%Branche ==      &
                  maille_r(imail-1)%Branche .and. &
                  abs(maille_r(imail)%AbscisseDeb-maille_r(imail - 1)%AbscisseFin).gt.EPS15) then
                  Erreur%Numero = 375
                  Erreur%ft     = err_375
                  Erreur%ft_c   = err_375c
                  call TRAITER_ERREUR( Erreur , 'du maillage par series' , imail - 1 , imail )
                  return
               endif
            end do

         case( TYPE_MAILLAGE_SERIE_PROFIL )

            pathNode = 'parametresPlanimetrageMaillage/maillage/maillageClavier/nbPlages'
            line = xcasReader(unitNum, pathNode)
            read(unit=line, fmt=*) nb_maille

            if( nb_maille <= 0 ) then
               Erreur%Numero = 306
               Erreur%ft     = err_306
               Erreur%ft_c   = err_306c
               call TRAITER_ERREUR( Erreur , 'Nombre de zones du maillage' )
               return
            end if

            if(.not.associated(maille_e)) allocate( maille_e(nb_maille) , STAT = retour )
            if( retour /= 0 ) then
               Erreur%Numero = 5
               Erreur%ft     = err_5
               Erreur%ft_c   = err_5c
               call TRAITER_ERREUR( Erreur , 'maille_e' )
               return
            end if

            if(.not.associated(maille_r)) allocate( maille_r(0) , STAT = retour )
            if( retour /= 0 ) then
               Erreur%Numero = 5
               Erreur%ft     = err_5
               Erreur%ft_c   = err_5c
               call TRAITER_ERREUR( Erreur , 'maille_r' )
               return
            end if

            pathNode = 'parametresPlanimetrageMaillage/maillage/maillageClavier/num1erProfPlage'
            line = xcasReader(unitNum, pathNode)
            read(unit=line, fmt=*) (maille_e(k)%ProfilDeb, k=1,nb_maille)

            pathNode = 'parametresPlanimetrageMaillage/maillage/maillageClavier/numDerProfPlage'
            line = xcasReader(unitNum, pathNode)
            read(unit=line, fmt=*) (maille_e(k)%ProfilFin, k=1,nb_maille)

            pathNode = 'parametresPlanimetrageMaillage/maillage/maillageClavier/pasEspacePlage'
            line = xcasReader(unitNum, pathNode)
            read(unit=line, fmt=*) (maille_e(k)%Pas, k=1,nb_maille)

            RCHID_2: DO K = 1, SIZE(PROFDEBBIEF)
               DO IMAIL = 1, NB_MAILLE
                  IF (MAILLE_E(IMAIL)%PROFILDEB .EQ. PROFDEBBIEF(K)) GOTO 300
               ENDDO
               ERREUR%NUMERO = 801
               ERREUR%FT = ERR_801
               ERREUR%FT_C = ERR_801C
               CALL TRAITER_ERREUR(ERREUR, 'NUMERO DU PREMIER PROFIL DE LA SERIE', K)
               RETURN
300            DO IMAIL = 1, NB_MAILLE
                  IF (MAILLE_E(IMAIL)%PROFILFIN .EQ. PROFFINBIEF(K)) CYCLE RCHID_2
               ENDDO
               ERREUR%NUMERO = 801
               ERREUR%FT = ERR_801
               ERREUR%FT_C = ERR_801C
               CALL TRAITER_ERREUR(ERREUR, 'NUMERO DU DERNIER PROFIL DE LA SERIE', K)
               RETURN
            ENDDO RCHID_2

            K = 1
            DO IMAIL = 1, NB_MAILLE

               if((MAILLE_E(IMAIL)%PROFILDEB .LT. 1) .OR. (MAILLE_E(IMAIL)%PROFILDEB .GT. SIZE(PROFIL))) then
                  Erreur%Numero = 311
                  Erreur%ft     = err_311
                  Erreur%ft_c   = err_311c
                  call TRAITER_ERREUR( Erreur , 'du numero du profil de debut d''une zone du maillage' , IMAIL )
                  return
               end if

               if((MAILLE_E(IMAIL)%PROFILFIN .LT. 1) .OR. (MAILLE_E(IMAIL)%PROFILFIN .GT. SIZE(PROFIL))) then
                  Erreur%Numero = 311
                  Erreur%ft     = err_311
                  Erreur%ft_c   = err_311c
                  call TRAITER_ERREUR( Erreur , 'du numero du profil de fin d''une zone du maillage' , IMAIL )
                  return
               end if

               if( maille_e(IMAIL)%ProfilFin <= maille_e(IMAIL)%ProfilDeb ) then
                  Erreur%Numero = 314
                  Erreur%ft     = err_314
                  Erreur%ft_c   = err_314c
                  call TRAITER_ERREUR( Erreur , 'des zones de maillage' , IMAIL )
                  return
               end if

               if( maille_e(IMAIL)%Pas <= 0. ) then
                  Erreur%Numero = 374
                  Erreur%ft   = err_374
                  Erreur%ft_c = err_374c
                  call TRAITER_ERREUR( Erreur , 'de la valeur du pas  d''une zone de maillage' , IMAIL )
                  return
               end if

               IF (IMAIL .GT. 1) THEN
                  IF (MAILLE_E(IMAIL - 1)%PROFILFIN .EQ. PROFFINBIEF(K)) THEN
                     K = K + 1
                     IF (MAILLE_E(IMAIL)%PROFILDEB .NE. (MAILLE_E(IMAIL - 1)%PROFILFIN + 1)) THEN
                        ERREUR%NUMERO = 804
                        ERREUR%FT = ERR_804
                        ERREUR%FT_C = ERR_804C
                        CALL TRAITER_ERREUR(ERREUR, 'NUMERO DU PREMIER PROFIL DE LA SERIE', (IMAIL - 1), IMAIL)
                        RETURN
                     ENDIF
                  ELSE
                     IF (MAILLE_E(IMAIL)%PROFILDEB .NE. MAILLE_E(IMAIL - 1)%PROFILFIN) THEN
                        ERREUR%NUMERO = 804
                        ERREUR%FT = ERR_804
                        ERREUR%FT_C = ERR_804C
                        CALL TRAITER_ERREUR(ERREUR, 'NUMERO DU PREMIER PROFIL DE LA SERIE', (IMAIL - 1), IMAIL)
                        RETURN
                     ENDIF
                  ENDIF
               ENDIF

            end do

         case( TYPE_MAILLAGE_INDIVIDUELLE )
            if(.not.associated(maille_r)) allocate( maille_r(0) , STAT = retour )
            if( retour /= 0 ) then
               Erreur%Numero = 5
               Erreur%ft     = err_5
               Erreur%ft_c   = err_5c
               call TRAITER_ERREUR( Erreur , 'maille_e' )
               return
            end if

            if(.not.associated(maille_e)) allocate( maille_e(0) , STAT = retour )
            if( retour /= 0 ) then
               Erreur%Numero = 5
               Erreur%ft     = err_5
               Erreur%ft_c   = err_5c
               call TRAITER_ERREUR( Erreur , 'maille_e' )
               return
            end if

            pathNode = 'parametresPlanimetrageMaillage/maillage/maillageClavier/nbSections'
            line = xcasReader(unitNum, pathNode)
            read(unit=line, fmt=*) nb_section

            if( nb_section <= 0 ) then
               Erreur%Numero = 306
               Erreur%ft     = err_306
               Erreur%ft_c   = err_306c
               call TRAITER_ERREUR( Erreur , 'Nombre de sections de calcul du maillage' )
               return
            end if

            allocate( itab(nb_section) , STAT = retour )
            if( retour /= 0 ) then
               Erreur%Numero = 5
               Erreur%ft     = err_5
               Erreur%ft_c   = err_5c
               call TRAITER_ERREUR( Erreur , 'itab' )
               return
            end if

            allocate( rtab(nb_section) , STAT = retour )
            if( retour /= 0 ) then
               Erreur%Numero = 5
               Erreur%ft     = err_5
               Erreur%ft_c   = err_5c
               call TRAITER_ERREUR( Erreur , 'rtab' )
               return
            end if

            allocate( X(nb_section) , STAT = retour )
            if( retour /= 0 ) then
               Erreur%Numero = 5
               Erreur%ft     = err_5
               Erreur%ft_c   = err_5c
               call TRAITER_ERREUR( Erreur , 'X' )
               return
            end if

            pathNode = 'parametresPlanimetrageMaillage/maillage/maillageClavier/branchesSection'
            line = xcasReader(unitNum, pathNode)
            read(unit=line, fmt=*) itab

            section%Branche = itab(1)
            if(SECTION%BRANCHE .NE. 1) then
                Erreur%Numero = 325
                Erreur%ft     = err_325
                Erreur%ft_c   = err_325c
                call TRAITER_ERREUR( Erreur , 1)
                return
            end if

            pathNode = 'parametresPlanimetrageMaillage/maillage/maillageClavier/absSection'
            line = xcasReader(unitNum, pathNode)
            read(unit=line, fmt=*) rtab

            RCHID_3: DO K = 1, SIZE(PROFDEBBIEF)
               DO IMAIL = 1, NB_SECTION
                  IF (ITAB(IMAIL) .EQ. K) GOTO 400
               ENDDO
               ERREUR%NUMERO = 800
               ERREUR%FT = ERR_800
               ERREUR%FT_C = ERR_800C
               CALL TRAITER_ERREUR(ERREUR, 'BRANCHES DES SECTIONS DE CALCUL', K)
               RETURN
400            DO IMAIL = 1, NB_SECTION
                  IF ((ITAB(IMAIL) .EQ. K) .AND. (RTAB(IMAIL) .EQ. ABSCRELEXTDEBBIEF(K))) GOTO 500
               ENDDO
               ERREUR%NUMERO = 801
               ERREUR%FT = ERR_801
               ERREUR%FT_C = ERR_801C
               CALL TRAITER_ERREUR(ERREUR, 'ABSCISSES DES SECTIONS DE CALCUL', K)
               RETURN
500            DO IMAIL = 1, NB_SECTION
                  IF ((ITAB(IMAIL) .EQ. K) .AND. (RTAB(IMAIL) .EQ. ABSCRELEXTFINBIEF(K))) CYCLE RCHID_3
               ENDDO
               ERREUR%NUMERO = 801
               ERREUR%FT = ERR_801
               ERREUR%FT_C = ERR_801C
               CALL TRAITER_ERREUR(ERREUR, 'ABSCISSES DES SECTIONS DE CALCUL', K)
               RETURN
            ENDDO RCHID_3

            section%AbscisseRel = rtab(1)
            if( section%AbscisseRel < AbscRelExtDebBief(section%Branche) .or. &
                section%AbscisseRel > AbscRelExtFinBief(section%Branche)) then
                Erreur%Numero   = 334
                Erreur%ft       = err_334
                Erreur%ft_c     = err_334c
                call TRAITER_ERREUR( Erreur , 'fin' , 1, SECTION%BRANCHE)
                return
           endif

            ! Passage en abscisses absolues
            X(1) = ABS_ABS_S            ( &
                section%Branche         , &
                section%AbscisseRel     , &
                Profil                  , &
                ProfDebBief             , &
                ProfFinBief             , &
                Erreur                    &
                                      )
            if( Erreur%Numero /= 0 ) then
               return
            end if

            do k = 2 , nb_section

               section%Branche     = itab(k)
               section%AbscisseRel = rtab(k)

               IF ((SECTION%BRANCHE .LT. 1) .OR. (SECTION%BRANCHE .GT. SIZE(PROFDEBBIEF))) THEN
                  ERREUR%NUMERO = 802
                  ERREUR%FT = ERR_802
                  ERREUR%FT_C = ERR_802C
                  CALL TRAITER_ERREUR(ERREUR, 'BRANCHES DES SECTIONS DE CALCUL', K)
                  RETURN
               ENDIF
               IF ((SECTION%ABSCISSEREL .LT. ABSCRELEXTDEBBIEF(SECTION%BRANCHE)) .OR. &
               (SECTION%ABSCISSEREL .GT. ABSCRELEXTFINBIEF(SECTION%BRANCHE))) THEN
                  ERREUR%NUMERO = 334
                  ERREUR%FT = ERR_334
                  ERREUR%FT_C = ERR_334C
                  CALL TRAITER_ERREUR(ERREUR, 'ABSCISSES DES SECTIONS DE CALCUL', K, SECTION%BRANCHE)
                  RETURN
               ENDIF

               ! Passage en abscisses absolues
               X(k) = ABS_ABS_S            ( &
                   section%Branche         , &
                   section%AbscisseRel     , &
                   Profil                  , &
                   ProfDebBief             , &
                   ProfFinBief             , &
                   Erreur                    &
                                      )
               if( Erreur%Numero /= 0 ) then
                  return
               end if

               if( X(k) <= X(k-1) ) then
                  Erreur%Numero = 310
                  Erreur%ft     = err_310
                  Erreur%ft_c   = err_310c
                  call TRAITER_ERREUR( Erreur , k )
                  return
               end if
            end do

            deallocate(itab)
            deallocate(rtab)

      case default

         allocate( maille_r(0) , STAT = retour )
         if( retour /= 0 ) then
            Erreur%Numero = 5
            Erreur%ft     = err_5
            Erreur%ft_c   = err_5c
            call TRAITER_ERREUR( Erreur , 'maille_e' )
            return
         end if

         allocate( maille_e(0) , STAT = retour )
         if( retour /= 0 ) then
            Erreur%Numero = 5
            Erreur%ft     = err_5
            Erreur%ft_c   = err_5c
            call TRAITER_ERREUR( Erreur , 'maille_e' )
            return
         end if

         !----------------------------------------------------
         ! cas TYPE_MAILLAGE_PRECEDENT et TYPE_MAILLAGE_PROFIL
         !----------------------------------------------------

      end select

   endif

   !-------------------
   ! Calcul du maillage
   !-------------------
   call MAILLER        ( &
        X              , & ! maillage
        maille_r       , & ! mailles reelles
        maille_e       , & ! mailles entieres
        Profil         , & ! tableau des profils
        TypeMaillage   , & ! type de calcul du maillage
        impression_geo , & ! test d'impression de la geometrie
        UniteListing   , & !
        Erreur           & ! Erreur
                       )
   if( Erreur%numero /= 0 ) then
      return
   endif

   ! Controle sur la taille du maillage
   !-----------------------------------
   !--------------------------------
   ! Sauvegarde eventuel du maillage
   !--------------------------------
   pathNode = 'parametresPlanimetrageMaillage/maillage/sauvMaillage'
   line = xcasReader(unitNum, pathNode)
   read(unit=line, fmt=*) sauvegarde_maillage

   if( sauvegarde_maillage ) then
      if (UniteListing >0) write(UniteListing,10490) 'OUI'
   else
      if (UniteListing >0) write(UniteListing,10490) 'NON'
   endif

   deallocate( maille_r , STAT = retour )
   if( retour /= 0 ) then
      Erreur%Numero = 6
      Erreur%ft   = err_6
      Erreur%ft_c = err_6c
      call TRAITER_ERREUR  (Erreur, 'maille_r')
      return
   end if

   deallocate (maille_e, STAT = retour)
   if( retour /= 0 ) then
      Erreur%Numero = 6
      Erreur%ft   = err_6
      Erreur%ft_c = err_6c
      call TRAITER_ERREUR  (Erreur, 'maille_e')
      return
   end if

   ! Fin des traitements
   !Erreur%arbredappel = !arbredappel_old

   return

   ! Formats d'ecriture

   10465 format (/,'MAILLAGE',/, &
               &  '---------',/)
   10470 format ('Type de calcul du maillage : ',A)
   10480 format ('Mode de saisie du maillage : ',A)
   10490 format (/,'Sauvegarde du maillage : ',A3)

   contains

   subroutine xerror(Erreur)

       use M_MESSAGE_C
       use M_ERREUR_T            ! Type ERREUR_T

       type(ERREUR_T)                   , intent(inout) :: Erreur

       Erreur%Numero = 704
       Erreur%ft     = err_704
       Erreur%ft_c   = err_704c
       call TRAITER_ERREUR( Erreur )

       return

   end subroutine xerror

end subroutine CALC_MAILLAGE
