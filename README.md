Hall-A Single Arm Monte Carlo Simulation Tool in C++
====================================================

History
-------

 * 11/16/2014 - Started updating code. Adding comments and making code more readable.
 * 9/30/2013  - Release by Zhihong Ye, but most credits are given to Alexander Deur and Huan Yao. 

About
-----

This is a C++ version of SAMC which was orginally developed by A. Deur using FORTRAN,
and was converted into a C++ package by Huan Yao. I made some small modification and added a generator. 
The package is released "AS IT IS" and we (especially A. Deur and H. Yao) don't hold any responsibility 
for updating the code, removing bugs, etc. Make sure you know the code well and feel free to do any 
modifications. 

It is a well-organized package and you can use it by only changing the input file, where all parameters
related to the experimental setup is defined (e.g., c12_input.dat). The generated events are stored in 
a ROOT file. The event seeds can be uniformly generated in the code, or can be specified outside the code 
using a generator (no physics yet) in ./Generator if you have special request.

The HRS transportation functions were generated using SNACK by J. LeRose and they are coded in the 
FORTRAN subroutines.


