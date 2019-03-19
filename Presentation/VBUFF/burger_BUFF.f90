!Mercier Wilfried
!09/01/2019
!Equation de burger - VBUFF
!Version 7 avec des tableaux de vitesse privés et buffers entre

PROGRAM Burger
   USE OMP_LIB
   IMPLICIT NONE
   
   !...Declaration...
      !...Variables...
      REAL(KIND=8)                                 :: stept                !Pas selon t
      REAL(KIND=8)                                 :: stepx                !Pas selon x
      REAL(KIND=8)                                 :: fder                 !Derivee premiere de v /x
      REAL(KIND=8)                                 :: sder                 !Derivee seconde de v /x
      INTEGER                                      :: x, t                 !Temps et espace
      INTEGER                                      :: NSTEPT               !Nombre de pas selon t
      INTEGER                                      :: NSSTEPX              !Nombre de pas dans les sous-tableaux 
      INTEGER                                      :: NSSSTEPX             !Nombre de pas dans les sous-tableaux
      INTEGER                                      :: NSSTEMP              !NSSTEPX+1
      INTEGER                                      :: NBUFF                !Nombre de buffers
      INTEGER                                      :: proc, nbprocs        !numero thread et nombre de threads
      integer                                      :: maxproc              !Numero thread max
      INTEGER                                      :: unitOut
      INTEGER                                      :: temp, tmp            !Stockage temporaire
      CHARACTER(LEN=30)                            :: filename
      LOGICAL                                      :: isOdd
      !...Parametres...
      REAL(KIND=8),  PARAMETER                     :: PI       = 4.0*atan(1.0_8) !PI
      REAL(KIND=8),  PARAMETER                     :: RE       = 20.0_8          !Nb de Reynolds
      REAL(KIND=8),  PARAMETER                     :: x0       = 0.0_8           !Premier bord
      REAL(KIND=8),  PARAMETER                     :: x1       = 1.0_8           !Second bord
      REAL(KIND=8),  PARAMETER                     :: TotTime  = 1.0_8           !Temps total
      INTEGER,       PARAMETER                     :: NSTEPX   = 1000            !Nombre de pas selon x
      INTEGER,       PARAMETER                     :: BUFFCLL  = 2*1             !Nombre de buffer par tableau
      !...Tableaux...
      REAL(KIND=8),  DIMENSION(:),  ALLOCATABLE    :: tab_x                !Tableau des positions
      REAL(KIND=8),  DIMENSION(:),  ALLOCATABLE    :: tab_v                !Tableau des vitesses
      REAL(KIND=8),  DIMENSION(:),  ALLOCATABLE    :: tab_temp             !Tableau temporaire
      REAL(KIND=8),  DIMENSION(:),  ALLOCATABLE    :: tab_buff             !Tableau buffer
      !...Fonctions...
      REAL(KIND=8)                                 :: OneDeriv             !Derivee premiere
      REAL(KIND=8)                                 :: TwoDeriv             !Derivee seconde
   !----------------------------------------------------------------------------------------------------------------------

   !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(tab_buff)
   
   !Initialisation pas de temps et espace
   stepx       =  (X1-X0)/(NSTEPX-1)
   stept       =  stepx*stepx

   !Nombre de processeurs
   nbprocs     =  OMP_GET_NUM_THREADS()
   maxproc     =  nbprocs-1

   !Dimensions des differents tableaux
   NSTEPT      =  INT(TotTime/stept)
   NSSTEPX     =  INT(NSTEPX/nbprocs)
   NSSSTEPX    =  INT(NSSTEPX+2)
   NSSTEMP     =  INT(NSSTEPX-1)
   NBUFF       =  INT(BUFFCLL*nbprocs)

   !Processeur fonctionnel
   proc        =  OMP_GET_THREAD_NUM()
   WRITE(filename, '(A8,I2.2)')"data/OMP", proc
   
   !Allocation tableaux vitesse et position
   ALLOCATE(tab_v(NSSTEPX))
   ALLOCATE(tab_temp(NSSTEPX))
   ALLOCATE(tab_x(NSSTEPX))

   !Allocation tableau buffer
   IF (proc .eq. 0) THEN
      ALLOCATE(tab_buff(NBUFF))
   ENDIF
   
   !Initialisation des positions
   DO x=1, NSSTEPX
      tab_x(x)  = x0 + ( proc*NSSTEPX + x-1)*stepx
   END DO

   !Initialisation des vitesses      
   DO x=1, NSSTEPX
      tab_v(x)    = SIN(2.0_8*PI*tab_x(x))
      tab_temp(x) = tab_v(x)
   END DO
      
   !Initialisation des buffers
   temp  =  proc*BUFFCLL
   DO x=1, BUFFCLL
      IF ( x .le. FLOOR(BUFFCLL/2.0_8) ) THEN 
         tab_buff(temp+x)  =  tab_v(x)
      ELSE
         tab_buff(temp+x)  =  tab_v(NSSTEPX-BUFFCLL+x)
      ENDIF
   END DO

   !Ecriture conditions initiales
   unitOut  =  10+proc
   OPEN(UNIT=unitOut, FILE=filename, STATUS='REPLACE')
   CALL Write2TabPlusToFile(tab_x, tab_v, proc, NSSTEPX, unitOut)
 
   !Calcul de l'evolution temporelle
   DO t=1, NSTEPT
   
      IF ( MOD(t, 2) .ne. 0 ) THEN
         isOdd =  .TRUE.
      ELSE
         isOdd =  .FALSE.
      ENDIF  

      !Mise a jour premiere case grace au buffer
      IF ( proc .eq. 0 ) THEN
         tmp   =  1
      ELSE
         tmp   =  proc*BUFFCLL
      ENDIF
      
      IF (isOdd) THEN
         fder     = OneDeriv(tab_temp(2), tab_temp(1), stepx, 0)
         sder     = TwoDeriv(tab_temp(2), tab_buff(tmp), tab_temp(1), stepx)
         tab_v(1) = stept*(sder/RE - tab_temp(1)*fder) + tab_temp(1)
      ELSE
         fder        = OneDeriv(tab_v(2), tab_v(1), stepx, 1)
         sder        = TwoDeriv(tab_v(2), tab_buff(tmp), tab_v(1), stepx)
         tab_temp(1) = stept*(sder/RE - tab_v(1)*fder) + tab_v(1)
      ENDIF

      !Mise a jour cases intermediaires sans buffer
      DO x=2, NSSTEMP, 1
         IF (isOdd) THEN
            fder        =  OneDeriv(tab_temp(x+1), tab_temp(x), stepx, 1)
            sder        =  TwoDeriv(tab_temp(x+1), tab_temp(x-1), tab_temp(x), stepx)
            tab_v(x)    =  stept*(sder/RE - tab_temp(x)*fder) + tab_temp(x)
         ELSE
            fder        =  OneDeriv(tab_v(x+1), tab_v(x), stepx, 1)
            sder        =  TwoDeriv(tab_v(x+1), tab_v(x-1), tab_v(x), stepx)     
            tab_temp(x) =  stept*(sder/RE - tab_v(x)*fder) + tab_v(x)                                
         ENDIF
      END DO

      !Mise a jour derniere case grace au buffer
      IF ( proc .eq. maxproc ) THEN
         tmp   =  NBUFF
      ELSE
         tmp   =  proc*BUFFCLL+3
      ENDIF
      
      IF (isOdd) THEN
         fder              =  OneDeriv(tab_buff(tmp), tab_temp(NSSTEPX), stepx, 1)
         sder              =  TwoDeriv(tab_buff(tmp), tab_temp(NSSTEMP), tab_temp(NSSTEPX), stepx)
         tab_v(NSSTEPX)    =  stept*(sder/RE - tab_temp(NSSTEPX)*fder) + tab_temp(NSSTEPX)
      ELSE
         fder              =  OneDeriv(tab_buff(tmp), tab_v(NSSTEPX), stepx, 1)
         sder              =  TwoDeriv(tab_buff(tmp), tab_v(NSSTEMP), tab_v(NSSTEPX), stepx)
         tab_temp(NSSTEPX) =  stept*(sder/RE - tab_v(NSSTEPX)*fder) + tab_v(NSSTEPX)
      ENDIF

      !Mise a jour buffer
      IF (isOdd) THEN
         temp  =  proc*BUFFCLL
         DO x=1, BUFFCLL
            IF ( x .le. FLOOR(BUFFCLL/2.0_8) ) THEN
               tab_buff(temp+x)  =  tab_v(x)
            ELSE
               tab_buff(temp+x)  =  tab_v(NSSTEPX-BUFFCLL+x)
            ENDIF
         END DO
      ELSE
         DO x=1, BUFFCLL
            IF ( x .le. FLOOR(BUFFCLL/2.0_8) ) THEN
               tab_buff(temp+x)  =  tab_temp(x)
            ELSE
               tab_buff(temp+x)  =  tab_temp(NSSTEPX-BUFFCLL+x)
            ENDIF
         END DO
      ENDIF

      DO x=1, NSSTEPX
         IF ( isOdd ) THEN
!            WRITE(*,*)"tab_v", x, tab_v(x)
         ELSE
!            WRITE(*,*)"tab_temp", x, tab_temp(x)
         ENDIF
      END DO

      !Ecrit dans fichier
      IF ( MOD(t, INT(NSTEPT/10)) .eq. 0 ) THEN
         IF ( MOD(t, 2) .ne. 0 ) THEN
            CALL Write2TabPlusToFile(tab_x, tab_v, t, NSSTEPX, unitOut)
         ELSE
            CALL Write2TabPlusToFile(tab_x, tab_temp, t, NSSTEPX, unitOut)
         ENDIF
      ENDIF
   END DO
   CLOSE(unitOut)

   !Desallocation tableaux
   DEALLOCATE(tab_v)
   DEALLOCATE(tab_temp)
   DEALLOCATE(tab_x)

   IF (proc .eq. 0) THEN
      DEALLOCATE(tab_buff)
   ENDIF

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
      WRITE(10, '(F13.7, 4X, F13.7)')tab1(i), tab2(i)
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
      WRITE(unitOut, '(I15, 4X, F13.7, 4X, F20.7)')time, tab1(i), tab2(i)
   END DO

END SUBROUTINE Write2TabPlusToFile
