The content of the folder
--------------------------
The folder contains MATLAB codes for the semi-global propagator, with examples,
the source codes of the paper from 2017, and some additional theoretical details on
the computation of the f_m(z,t) functions.

SemiGlobal.m implements the propagator in the general case.

SemiGlobalArnoldi_xp.m implements the propagator for quantum problems, with a Hamiltonian
which can be split into x dependent term and p dependent term. The Arnoldi approach is
employed for the computation of the function of matrix. This version was used to
obtain the results presented in the paper from 2017.

The other functions in the main folder are called from the main programs mentioned
above.

f_error.pdf contains some additional theoretical explanation on the computation of the
f_m(z,t) functions.


The content of the subfolder article_source_files
--------------------------------------------------
This folder contains the programs which compute the results presented in the paper, and 
additional data files.
The main folder has to be added to the MATLAB path in order that the program
SemiGlobalArnoldi_xp.m will be recognized by the calling functions in this folder.

The main program of the computation of the error decay curves of the semi-global propagator
is getSGdata4.m .

The main program of the computation of the error decay curves of Runge-Kutta of the 4'th 
order is RKerror.m .

Uexact_article.m computes the reference exact solution. The obtained solution is stored 
in the data file Uex_article.mat .

coulomb_optV240.mat is a data file which contains the required unperturbed potential
(Vabs240), the grid (x240), kinetic energy vector in the p domain (K240), the ground state
(fi0240), and the dipole moment which decays to 0 at the absorbing boundaries (xabs240). It
contains some additional data.

Vabs.mat is a data file which contains the optimized complex absorbing potential. This 
potential is added to the boundaries of the variable V0240 from coulomb_optV240.mat in
order to obtain the variable Vabs240 from coulomb_optV240.mat.


The content of the subfolder examples
--------------------------------------
The folder contains several examples for the application of SemiGlobal.m .
The main folder should be added to the MATLAB path in order that the program
SemiGlobal.m will be recognized by the calling functions in this folder.

test_harmonic.m tests the propagator for a forced harmonic oscillator problem.

testBECsg.m is similar to test_harmonic.m, but with the addition of a nonlinear BEC trap 
potential.

test_source_term.m is similar to test_harmonic.m, but with the addition of an arbitrary 
source term.

The subfolder article_results demonstrates the application of SemiGlobal.m
for obtaining the error decay curves presented in the paper. The main
program of the computation of the error decay curves is getSGdataSGcode.m .
The subfolder article_source_files has to be added to the MATLAB path in
order that the required data files will be recognized.
