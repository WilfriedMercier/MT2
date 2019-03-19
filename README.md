# MT2

## Overall description

Repository containing my 2nd Master project on parallel computing. The objective was to parrallellize a code which solved Burger equation on turbulent motion with both OpenMP and MPI. All the codes are in Fortran.

## Folders description

Folder | Description
:---: | :---
Scalaire | It contains the scalar version of the code which solved Burger equation
Vi (i in [3,8]) | Folders containing the different parrallellized versions of the code. To see how well the different versions behave, look at the SpeedUps.pdf file in Presentation folder
VMPI | The only MPI version avaiblable. The code can be compiled and launched but the results are not physically correct
Presentation | All the files related to the oral presentation are present here
