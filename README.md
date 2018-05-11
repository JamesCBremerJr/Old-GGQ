This package contains code for constructing generalized Gaussian quadrature
rules.  The following files comprise the package:

1.  The file print.f contains basic output routines.

2.  The file utils.f contains some simple utlity routines for sorting
and the like.

3.  The file qrsolve.f contains simple subroutines for solving systems of linear
equations.

4.  The file orthom.f contains simple subroutines for solving systems of
linear equations.

5.  The file gspiv.f contains code for conducting Gram-Schmidt orthogonalization


6.  The file legendre.f contains code for constructing Gauss-Legendre
quadratures and for manipulating Legendre expansions.

7.  The file adaptri.f contains code for adaptively integrating a user-supplied
function given over a user-specified triangle.

8.  The file legedisc.f contains code for representing functions as piecewise
Legendre expansions.

9.  The file newtls.f contains code for solving least squares problems.

10.  The file newton1d.f contains a one-dimensional ``generalized Gaussian''
quadrature code and the file newton2d.f contains the two-dimensional version
of the same.  Both codes rely on the utility routines in newtls.f.

11.  The file logquads.f contains code for constructing quadrature rules
for evaluating integrals of the form 

       1
    int   ( log|x-y| p(y) + q(y) ) dy
      -1

with p and q polynomials of specified degrees.