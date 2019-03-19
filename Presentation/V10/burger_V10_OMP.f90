!Mercier Wilfried
!09/01/2019
!Equation de burger - V10
!Ecriture dans plusieurs fichiers de maniere simultaneeburger_V7_OMP.f90

PROGRAM Burger
   USE OMP_LIB
   IMPLICIT NONE
   
   !...Declaration...
      !...Variables...
      REAL(KIND=8)                        :: stept                !Pas selon t
      REAL(KIND=8)                        :: stepx                !Pas selon x
      REAL(KIND=8)                        :: fder                 !Derivee premiere de v /x
      REAL(KIND=8)                        :: sder                 !Derivee seconde de v /x
      INTEGER                             :: x, t                 !Temps et espace
      INTEGER                             :: NSTEPT               !Nombre de pas selon t 
      INTEGER                             :: thread_num           !Numero du thread
      CHARACTER(LEN=50)                   :: fname                !Output filename
      !...Parametres...
      REAL(KIND=8), PARAMETER             :: PI = 4.0*atan(1.0_8) !PI
      REAL(KIND=8), PARAMETER             :: RE = 150.0_8         !Nb de Reynolds
      REAL(KIND=8), PARAMETER             :: x0 = 0.0_8           !Premier bord
      REAL(KIND=8), PARAMETER             :: x1 = 1.0_8           !Second bord
      REAL(KIND=8), PARAMETER             :: TotTime = 1.0_8      !Temps total
      INTEGER,      PARAMETER             :: NSTEPX = 2000        !Nombre de pas selon x
      !...Tableaux...
      REAL(KIND=8), DIMENSION(NSTEPX)     :: tab_x                !Tableau des positions
      REAL(KIND=8), DIMENSION(NSTEPX+2)   :: tab_v                !Tableau des vitesses
      REAL(KIND=8), DIMENSION(NSTEPX+2)   :: tab_temp             !Tableau temporaire
      !...Fonctions...
      REAL(KIND=8)                        :: OneDeriv             !Derivee premiere
      REAL(KIND=8)                        :: TwoDeriv             !Derivee seconde
   !-----------------------------------------------------------------

   !Initialisation pas de temps et espace
   stepx = (X1-X0)/(NSTEPX-1)
   stept = stepx*stepx
   NSTEPT   = INT(TotTime/stept)     

   !Initialisation des positions et des vitesses
   tab_v    = 0
   tab_temp = 0
   DO x=1, NSTEPX
      tab_x(x)       = x0 + (x-1)*stepx
      tab_v(x+1)     = SIN(2.0_8*PI*tab_x(x))
      tab_temp(x+1)  = tab_v(x+1)   
   END DO

   !Calcul de l'evolution temporelle
   !$OMP PARALLEL PRIVATE(fder, sder, fname, thread_num)
   thread_num  =  OMP_GET_THREAD_NUM()
   WRITE(fname, '(A8, I2.2)')"data/OMP", thread_num
   OPEN(UNIT=thread_num, FILE=fname, STATUS='REPLACE')
   CALL Write2TabPlusToFile(tab_x, tab_v, 0, NSTEPX, thread_num)
   DO t=1, NSTEPT
 
   !$OMP DO SCHEDULE(STATIC)
      DO x=2, NSTEPX+1
         !Calcul premiere et seconde derivees et met a jour le tableau
         IF ( MOD(t, 2) .ne. 0 ) THEN
            fder     = OneDeriv(tab_temp(x+1), tab_temp(x), stepx, 1)
            sder     = TwoDeriv(tab_temp(x+1), tab_temp(x-1), tab_temp(x), stepx)
            tab_v(x) = stept*(sder/RE - tab_temp(x)*fder) + tab_temp(x)
         ELSE
            fder        = OneDeriv(tab_v(x+1), tab_v(x), stepx, 1)
            sder        = TwoDeriv(tab_v(x+1), tab_v(x-1), tab_v(x), stepx)    
            tab_temp(x) = stept*(sder/RE - tab_v(x)*fder) + tab_v(x)                                         
         ENDIF
      END DO
   !$OMP END DO

   !$OMP MASTER
      !Ecrit dans fichier
      IF ( MOD(t, INT(NSTEPT/10)) .eq. 0 ) THEN
         IF ( MOD(t, 2) .ne. 0 ) THEN
            CALL Write2TabPlusToFile(tab_x, tab_v, t, NSTEPX, thread_num)
         ELSE
            CALL Write2TabPlusToFile(tab_x, tab_temp, t, NSTEPX, thread_num)
         ENDIF
      ENDIF
   !$OMP END MASTER
   END DO
   CLOSE(thread_num)
   !$OMP END PARALLEL

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
