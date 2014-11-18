Single Arm Monte Carlo 2
========================

Improvements
------------

Samc2 is a major rewrite and all the C++ interface has changed almost entirely. The original C++
code had many bugs, some more severe than others. 
Initially started as an attempt to make the code more readable with some bug fixes, it quickly evolved into a major rewrite.
This new version uses proper C++ techniques. Some of the major differences and changes include,

 * A good object model (as opposed to one object which does everthing)
 * Event class is actually a streamed event class for analysis. 
 * Readable code 
 * New build system using CMake. 
 * Directory structure separates source implementation and headers.
 * samc-config for build and linking.
 * ROOT libraries for scripting in addition to stand alone simulation
 * 

*NOTE*: This version is almost entirely different from the original SAMC. Please follow this documentation for 
    proper use

History
-------
 
 * 11/17/2014 - Original version had Undefined behavior! ... Started rewrite.
 * 11/16/2014 - Started updating code. Adding comments and making code more readable.
 * 9/30/2013  - Release by Zhihong Ye, but most credits are given to Alexander Deur and Huan Yao. 


About
-----

This is a *major* rewrite of the original SAMC fortran/c++ codes.

The HRS transportation functions were generated using SNACK by J. LeRose and they are coded in the 
FORTRAN subroutines.

Authors and Contributors
------------------------

This major rewrite
 * Whitney R. Armstrong (whit@temple.edu)    

Previous versions:
 * A. Deur (Original fortran version)
 * John LeRose  (SNAKE) HRS transport functions 
 * Huan Yao   (SAMC) First attempt to use C++.
 * Zhihong Ye (SAMC) Bug fixes 

