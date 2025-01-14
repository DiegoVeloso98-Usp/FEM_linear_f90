# FEM_Linear_f90
Linear FEM code, written in fortran90. 
2 Different solvers, SUPERLU and MA87

The code reads a 2D mesh in a specefic format.
 - Allows for linear, quadratic and cubic approximations.
 - It allows for particles generations inside the global mesh. The particles then contribute to the stiffness of the elements where they are placed. The particle generation is ramdom.
 - After the solution, it outputs 2 different files, one which can be viewed with AcadView (Visualization software from SET-EESC, University of São Paulo (USP) )
