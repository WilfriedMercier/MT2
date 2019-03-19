!Mercier Wilfried
!07/01/2019
!Equation de burger - V1
!Grille initialisee et fonction gerant la sortie ecrite

PROGRAM Burger
   IMPLICIT NONE
   
   !...Declaration...
      !...Variables...
      REAL(KIND=8)                     :: t, x                 !Temps et espace
      REAL(KIND=8)                     :: stept                !Pas selon t
      REAL(KIND=8)                     :: stepx                !Pas selon x
      INTEGER                          :: i                    !Variable d'incrementation
      !...Parametres...
      REAL(KIND=8), PARAMETER          :: PI = 4.0*atan(1.0_8) ! PI
      REAL(KIND=8), PARAMETER          :: RE = 1.0_8           !Nb de Reynolds
      REAL(KIND=8), PARAMETER          :: x0 = 0.0_8           !Premier bord
      REAL(KIND=8), PARAMETER          :: x1 = 1.0_8           !Second bord
      INTEGER,      PARAMETER          :: NSTEPX = 100         !Nombre de pas selon x
      !...Tableaux...
      REAL(KIND=8), DIMENSION(nstepx)  :: tab_x                !Tableau des positions
      REAL(KIND=8), DIMENSION(nstepx)  :: tab_v                !Tableau des vitesses
   !-----------------------------------------------------------------

   !Initialisation pas de temps et espace
   stepx = REAL( (X1-X0)/(NSTEPX-1) )
   stept = REAL( stepx*stepx )

   !Initialisation des positions et des vitesses
   DO i=1, nstepx
      tab_x(i) = x0 + (i-1)*stepx
      tab_v(i) = SIN(2.0_8*PI*tab_x(i))
   END DO

   !Ecriture conditions initiales dans un fichier
   CALL WriteTabToFile(tab_x, tab_v, NSTEPX, "data_init")



END PROGRAM Burger



!----------------------------------------------------------------------------------------------------------------------
!>Ecrit dans un fichier deux tableaux 1D
SUBROUTINE WriteTabToFile(tab1, tab2, length, fname)
   IMPLICIT NONE

   !...Declaration...
      !...Variables locales...
      INTEGER,                    INTENT(IN) :: length   !Taille tableau
      REAL(KIND=8), DIMENSION(*), INTENT(IN) :: tab1     !Tableau
      REAL(KIND=8), DIMENSION(*), INTENT(IN) :: tab2     !Tableau
      CHARACTER(LEN=*),           INTENT(IN) :: fname
      !...non locales...
      INTEGER                                :: i
      
   !Ouverture fichier
   OPEN(UNIT=10, FILE=fname, STATUS='REPLACE')

   !Ecriture fichier
   DO i=1, length
      WRITE(10, '(F10.8, 4X, F10.8)')tab1(i), tab2(i)
   END DO

   !Fermeture fichier
   CLOSE(10)

END SUBROUTINE WriteTabToFile
