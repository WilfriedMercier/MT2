!Mercier Wilfried
!09/01/2019
!Equation de burger - V8
!Essaie optimisation du code

PROGRAM Burger
   IMPLICIT NONE
   
   !...Declaration...
      !...Variables...
      REAL(KIND=8)                     :: stept                   !Pas selon t
      REAL(KIND=8)                     :: stepx                   !Pas selon x
      REAL(KIND=8)                     :: fder                    !Derivee premiere de v /x
      REAL(KIND=8)                     :: sder                    !Derivee seconde de v /x
      INTEGER                          :: x, t                    !Temps et espace
      INTEGER                          :: NSTEPT                  !Nombre de pas selon t 
      INTEGER                          :: temp                    !Stockage temporaire
      INTEGER                          :: step                    !Decalge pour parcourir le bon tableau
      LOGICAL                          :: IsEven                  !Variable logique pour savoir quel tableau utiliser
      !...Parametres...
      REAL(KIND=8), PARAMETER          :: PI = 4.0*atan(1.0_8)    !PI
      REAL(KIND=8), PARAMETER          :: RE = 150.0_8            !Nb de Reynolds
      REAL(KIND=8), PARAMETER          :: x0 = 0.0_8              !Premier bord
      REAL(KIND=8), PARAMETER          :: x1 = 1.0_8              !Second bord
      REAL(KIND=8), PARAMETER          :: TotTime  = 1.0_8        !Temps total
      INTEGER,      PARAMETER          :: NSTEPX   = 2000         !Nombre de pas selon x
      INTEGER,      PARAMETER          :: NTOT     = 2*NSTEPX+2   !Nombre de pas selon x
      !...Tableaux...
      REAL(KIND=8), DIMENSION(NTOT)    :: tab_x                   !Tableau des positions
      REAL(KIND=8), DIMENSION(NTOT)    :: tab_v                   !Tableau des vitesses
      !...Fonctions...
      REAL(KIND=8)                     :: OneDeriv                !Derivee premiere
      REAL(KIND=8)                     :: TwoDeriv                !Derivee seconde
   !-----------------------------------------------------------------

   !Initialisation pas de temps et espace
   stepx    = (X1-X0)/(NSTEPX-1)
   stept    = stepx*stepx
   NSTEPT   = INT(TotTime/stept)     

   !Initialisation des positions
   tab_x(1) =  0.0_8
   DO x=0, NSTEPX-1, 1
      temp           =  2*(x+1)
      tab_x(temp)    =  x0 + x*stepx
      tab_x(temp+1)  =  tab_x(temp)
   END DO
   tab_x(NTOT) =  0.0_8

   !Initialisation vitesses
   tab_v(1) = 0.0_8
   DO x=2, NTOT-1
      tab_v(x)    = SIN(2.0_8*PI*tab_x(x))
   END DO
   tab_v(NTOT)    = 0.0_8

   !Calcul de l'evolution temporelle
   OPEN(UNIT=11, FILE="data/OMP_datafile", STATUS='REPLACE')
   CALL WriteOutputToFile(tab_x, tab_v, 0, NTOT, 11, .TRUE. )

   DO t=1, NSTEPT
   
      IF ( MOD(t, 2) .ne. 0 ) THEN
         IsEven   = .FALSE.
         step     = 1
         
         fder     = OneDeriv(tab_v(4), tab_v(1), stepx, 1)
         sder     = TwoDeriv(tab_v(4), tab_v(1), tab_v(2), stepx)
         tab_v(3) = stept*(sder/RE - tab_v(2)*fder) + tab_v(2)
      ELSE
         IsEven   = .TRUE.
         step     = 0
         
         fder     = OneDeriv(tab_v(5), tab_v(1), stepx, 1)
         sder     = TwoDeriv(tab_v(5), tab_v(1), tab_v(3), stepx)
         tab_v(2) = stept*(sder/RE - tab_v(3)*fder) + tab_v(3)
      ENDIF

      DO x=4+step, NTOT-3, 2
         !Calcul premiere et seconde derivees
         IF (IsEven) THEN
            fder     = OneDeriv(tab_v(x+3), tab_v(x-1), stepx, 1)
            sder     = TwoDeriv(tab_v(x+3), tab_v(x-1), tab_v(x+1), stepx)
            tab_v(x) = stept*(sder/RE - tab_v(x+1)*fder) + tab_v(x+1)
         ELSE
            fder     = OneDeriv(tab_v(x+1), tab_v(x-3), stepx, 1)
            sder     = TwoDeriv(tab_v(x+1), tab_v(x-3), tab_v(x-1), stepx)
            tab_v(x) = stept*(sder/RE - tab_v(x-1)*fder) + tab_v(x-1)
         ENDIF
      END DO
      
      IF (IsEven) THEN
         temp        =  NTOT-3
         x           =  NTOT-1
         fder        =  OneDeriv(tab_v(NTOT), tab_v(temp), stepx, 1)
         sder        =  TwoDeriv(tab_v(NTOT), tab_v(temp), tab_v(x), stepx)
         tab_v(x-1)  =  stept*(sder/RE - tab_v(x)*fder) + tab_v(x)
      ELSE
         temp        =  NTOT-4
         x           =  NTOT-2
         fder        =  OneDeriv(tab_v(NTOT), tab_v(temp), stepx, 1)
         sder        =  TwoDeriv(tab_v(NTOT), tab_v(temp), tab_v(x), stepx)
         tab_v(x+1)  =  stept*(sder/RE - tab_v(x)*fder) + tab_v(x)
      ENDIF
      
      !Ecrit dans fichier
      IF ( MOD(t, INT(NSTEPT/10)) .eq. 0 ) THEN
         CALL WriteOutputToFile(tab_x, tab_v, t, NTOT, 11, IsEven)
      ENDIF
   END DO
   CLOSE(11)

END PROGRAM Burger


!---------------------------------------------------------------------------------------------------------------------
!>Derivee premiere non-centree/centree
REAL(KIND=8) FUNCTION OneDeriv(valSup, val, step, isCenteredFlag)
   IMPLICIT NONE

   !...non locales...
   REAL(KIND=8), INTENT(IN) :: valSup           !Valeur a i+1
   REAL(KIND=8), INTENT(IN) :: val              !Valeur a i
   REAL(KIND=8), INTENT(IN) :: step             !Pas
   INTEGER,      INTENT(IN) :: isCenteredFlag   !0 si non centre, 1 si centre

   OneDeriv = (valSup-val)/step
!   IF (OneDeriv .gt. HUGE(OneDeriv) ) THEN
!      WRITE(*,*)valSup, val, step
!   ENDIF
   
   IF (isCenteredFlag .EQ. 1) THEN
      OneDeriv = OneDeriv/2.0_8
   ENDIF
   
END FUNCTION OneDeriv


!---------------------------------------------------------------------------------------------------------------------
!>Derivee premiere ordre deux
REAL(KIND=8) FUNCTION OneDerivSecondOrder(valSupSup, valSup, val, step)
   IMPLICIT NONE

   !...non locales...
   REAL(KIND=8), INTENT(IN) :: valSupSup  !Valeur a i+2
   REAL(KIND=8), INTENT(IN) :: valSup     !Valeur a i+1
   REAL(KIND=8), INTENT(IN) :: val        !Valeur a i
   REAL(KIND=8), INTENT(IN) :: step       !Pas

   OneDerivSecondOrder = (4.0_8*valSup - valSupSup - 3.0_8*val)/(2.0_8*step)

END FUNCTION OneDerivSecondOrder


!---------------------------------------------------------------------------------------------------------------------
!>Derivee seconde
REAL(KIND=8) FUNCTION TwoDeriv(valSup, valMin, val, step)
   IMPLICIT NONE

   !...non locales...
   REAL(KIND=8), INTENT(IN) :: valSup  !Valeur a i+1
   REAL(KIND=8), INTENT(IN) :: valMin  !Valeur a i-1
   REAL(KIND=8), INTENT(IN) :: val     !Valeur a i
   REAL(KIND=8), INTENT(IN) :: step    !Pas

   TwoDeriv = REAL( (valSup+valMin-2.0_8*val)/(step*step)  )

END FUNCTION TwoDeriv


!----------------------------------------------------------------------------------------------------------------------
!>Ecrit dans un fichier deux tableaux 1D
SUBROUTINE Write2TabToFile(tab1, tab2, length, fname)
   IMPLICIT NONE

   !...Declaration...
      !...non locales...
      INTEGER,                    INTENT(IN) :: length   !Taille tableau
      REAL(KIND=8), DIMENSION(*), INTENT(IN) :: tab1     !Tableau
      REAL(KIND=8), DIMENSION(*), INTENT(IN) :: tab2     !Tableau
      CHARACTER(LEN=*),           INTENT(IN) :: fname
      !...Variables locales...
      INTEGER                                :: i
      
   !Ouverture fichier
   OPEN(UNIT=10, FILE=fname, STATUS='REPLACE')

   !Ecriture fichier
   DO i=1, length
      WRITE(10, '(F10.8, 4X, F10.8)')tab1(i), tab2(i)
   END DO

   !Fermeture fichier
   CLOSE(10)

END SUBROUTINE Write2TabToFile


!--------------------------------------------------------------------------------------------------------$
!>Ecrit dans un fichier deux tableaux 1D
SUBROUTINE Write2TabPlusToFile(tab1, tab2, time, length, unitOut)
   IMPLICIT NONE

   !...Declaration...
   !...non locales...
   INTEGER,                    INTENT(IN) :: length, unitOut, time   !Taille tableau et unite sortie
   REAL(KIND=8), DIMENSION(*), INTENT(IN) :: tab1, tab2              !Tableaux
   !...Variables locales...
   INTEGER                                :: i

   !Ecriture fichier
   DO i=1, length
      WRITE(unitOut, '(I15, 4X, F13.7, 4X, F13.7)')time, tab1(i), tab2(i)
   END DO

END SUBROUTINE Write2TabPlusToFile


!--------------------------------------------------------------------------------------------------------$
!>Ecrit dans un fichier deux tableaux 1D
SUBROUTINE WriteOutputToFile(tab1, tab2, time, length, unitOut, isEven)
   IMPLICIT NONE

   !...Declaration...
   !...non locales...
   INTEGER,                    INTENT(IN) :: length, unitOut, time   !Taille tableau et unite sortie
   LOGICAL,                    INTENT(IN) :: isEven
   REAL(KIND=8), DIMENSION(*), INTENT(IN) :: tab1, tab2              !Tableaux
   !...Variables locales...
   INTEGER                                :: i, start
   
   IF (isEven) THEN
      start =  2
   ELSE
      start =  3
   ENDIF

   !Ecriture fichier
   WRITE(unitOut, '(I15, 4X, F13.7, 4X, F13.7)')time, tab1(1), tab2(1)
   DO i=start, length-1, 2
      WRITE(unitOut, '(I15, 4X, F13.7, 4X, F13.7)')time, tab1(i), tab2(i)
   END DO
   WRITE(unitOut, '(I15, 4X, F13.7, 4X, F13.7)')time, tab1(length), tab2(length)

END SUBROUTINE WriteOutputToFile


