# A MATLAB interface for L-BFGS-B

###Update

As of March 24, 2014, the MATLAB code supports the latest version of
the L-BFGS-B solver (version 3.0), and is compatible with
[GNU Octave](https://www.gnu.org/software/octave). Thank you to
[Jose Vallet](http://github.com/josombio) for providing these updates.

###Overview

L-BFGS-B is a collection of Fortran 77 routines for solving
large-scale nonlinear optimization problems with bound constraints on
the variables. One of the key features of the nonlinear solver is that
knowledge of the Hessian is not required; the solver computes search
directions by keeping track of a quadratic model of the objective
function with a limited-memory BFGS (Broyden-Fletcher-Goldfarb-Shanno)
approximation to the Hessian *(see Note #1)*. The algorithm was
developed by Ciyou Zhu, Richard Byrd and Jorge Nocedal. For more
information, go to the
[original distribution site](http://www.ece.northwestern.edu/~nocedal/lbfgsb.html)
for the L-BFGS-B software package.

I've designed an interface to the L-BFGS-B solver so that it can be
called like any other function in MATLAB *(see Note #2)*. See the text
below for more information on installing and calling this function in
MATLAB. Along the way, I've also developed a C++ class that
encapsulates all the "messy" details in executing the L-BFGS-B
code. See below for instructions on how to incorporate this C++ class
info your own code.

This code has been tested in MATLAB and Octave on Linux, and is
integrated into [GPML](http://www.gaussianprocess.org/gpml/code). It
has also been executed successfully in Windows; see below.

###License

Copyright (c) 2013, 2014 Peter Carbonetto

The lbfgsb-matlab project repository is free software: you can
redistribute it and/or modify it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html) as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
**without any warranty**; without even the implied warranty of
**merchantability** or **fitness for a particular purpose**. See
[LICENSE](LICENSE) for more details.

###Installation

Follow these steps to compile and install **lbfgsb-matlab** for MATLAB
and Octave. Some of these steps are applicable to MATLAB or Octave
only, and you can skip them if you want to install the interface for
only one of them. These installation instructions assume you have
Linux, although similar steps should also work for other UNIX-based
operating systems, such as Mac OS X. These instructions also assume
you have [GNU Make](http://www.gnu.org/software/make) installed. (We
have successfully compiled this source code with gcc versions 4.6.3
and 4.6.4, under Ubuntu Linux 12.10 and 13.10, and have successfully
run the MATLAB scripts in Octave 3.6.2, 3.6.4, 3.8.0 and 3.8.1.)

**For Windows users:** If you have a Windows operating system,
Guillaume Jacquenot has generously provided
[instructions](INSTALL.WINDOWS.txt) for compiling the software on
Windows with [gnumex](http://gnumex.sourceforge.net) and Mingw. These
instructions are slightly out of date, and suggestions for improving
these installation instructions are welcome.

We will create a MEX file, which is basically a file that contains a
routine that can be called from MATLAB as if it were a built-it
function. To learn about MEX files, I refer you to
[this document](http://www.mathworks.com/support/tech-notes/1600/1605.html)
on the MathWorks website.

**Download the source code.** Clone or fork the this repository, or
download this repository as a
[ZIP file](http://github.com/pcarbo/lbfgsb-matlab/archive/master.zip)
and unpack the ZIP file. The [src](src) subdirectory contains some
MATLAB files (ending in .m), some C++ source and header files (ending
in .h and .cpp), a single Fortran 77 source file,
[solver.f](src/solver.f), containing the L-BFGS-B routines, and a
[Makefile](src/Makefile). I've included a version of the Fortran
routines that is more recent than what is available for download at
the
[distribution site](http://www.ece.northwestern.edu/~nocedal/lbfgsb.html)
at Northwestern University.

In the following, we will refer to the directory containing your local
copy of the repository as LBFGSB_HOME. Wherever you find LBFGSB_HOME,
you should substitute it with your own directory.

**Install the C++ and Fortran 77 compilers.** In order to build the
MEX file, you will need both a C++ compiler and a Fortran 77
compiler. Unfortunately, you can't just use any compiler. You have to
use the precise one supported by MATLAB. For instance, if you are
running Mac OS X 7.2 on a Linux operating system, you will need to
install [GNU compiler collection](http://www.gnu.org/software/gcc)
(GCC) version 3.4.4. Even if you already have a compiler installed on
Linux, it may be the wrong version, and if you use the wrong version
things could go horribly wrong. It is important that you use the
correct version of the compiler, otherwise you will encounter linking
errors. Refer to
[this document](http://www.mathworks.com/support/tech-notes/1600/1601.html)
to find out which compiler is supported by your version of MATLAB.

You may also be able to use the **apt-get** program to install
gfortran, as in:

    sudo apt-get install build-essential gfortran

**Install Octave.** You will of course need a working installation of
Matlab or Octave. In the case of Octave, if you don't have it the
easiest way to install it is from the repositories using **apt-get**:

    sudo apt-get install octave

If you want to use a more recent version of Octave, you will have to
[install Octave from the source code](http://www.gnu.org/software/octave/doc/interpreter/Installation.html).
After that, find the path to the mkoctfile compiler (in the bin
directory) and the include directory. You will need those paths for
steps below.

Additionally you will need the Octave development tools, libs and
header files to compile MEX files:

    sudo apt-get install liboctave-dev

**Configure MATLAB.** Next, you need to set up and configure MATLAB
to build MEX Files. This is explained quite nicely in
[this MathWorks support document](http://www.mathworks.com/support/tech-notes/1600/1605.html).

**Simple compilation procedure for MATLAB.** (Thanks to Stefan
Harmeling.) You might be able to build and compile the MEX File from
source code simply by calling MATLAB's mex program with the C++ and
Fortran source files as inputs, like so:

    mex -output lbfgsb *.cpp solver.f

If that doesn't work (e.g. you get linking errors, or it makes MATLAB
crash when you try to call it), then you may have to follow the more
complicated installation instructions below.

**Modify the Makefile.** You are almost ready to build the MEX
file. But before you do so, you need to edit the Makefile to coincide
with your system setup. Edit the Makefile located inside the
LBFGSB_HOME/src directory, following the instructions provided in this
file. With this Makefile you may create:

- A mex file for Matlab using the mex compiler (mex taget).
- A mex file for Matlab using your C compiler (nomex target).
- A mex file for Octave using mkoctfile (oct target).

Regarding the compilation for Octave, the variable **OCTAVE_INCLUDE**
is the directory where the Octave header files required for
development are installed. Make sure that it points to the correct
directory corresponding to your Octave installation. **OCT** is the
command to call the mkoctfile compiler. If mkoctfile is in the path,
simply leave it as it is. If not, include in the variable the full
path to it. If you have several versions of Octave installed, make
sure that OCT points to the mkoctfile program that corresponds to the
Octave version you would like to use. 

For MATLAB, the variable **MEX** is the executable called
to build the MEX file. The variables **CXX** and **F77** must be the
names of your C++ and Fortran 77 compilers. **MEXSUFFIX** is the MEX
file extension specific to your platform. (For instance, the extension
for Linux is **mexglx**.) The variable **MATLAB_HOME** must be the
base directory of your MATLAB installation. Finally, **CFLAGS** and
**FFLAGS** are options passed to the C++ and Fortran compilers,
respectively. Often these flags coincide with your MEX options file
(see
[here](http://www.mathworks.com/support/tech-notes/1600/1605.html) for
more information on the options file). But often the settings in the
MEX option file are incorrect, and so the options must be set
manually. MATLAB requires some special compilation flags for various
reasons, one being that it requires position-independent code. These
instructions are vague (I apologize), and this step may require a bit
of trial and error until before you get it right.

**INSTALLDIR** is the installation directory relative to the
LBFGSB_HOME/src directory (see Installation below).

It may be helpful to look at the GCC documentation in order to
understand what these various compiler flags mean.

**Build the MEX file.** If you are in the directory containing all the
source files, calling **make** at the command prompt will first
compile the Fortran and C++ source files into object code (.o
files). After that, the make program calls the MEX script, which links
all the object files together into a single MEX file. Call **make**
from the LBFGSB_HOME/src directory using the appropriate target, one
of:

    make oct
    make nomex 
    make mex

If you didn't get any errors, then you are ready to try out the
bound-constrained solver in MATLAB or Octave. Note that even if you
didn't get any errors, there's still a possibility that you didn't
link the MEX file properly, in which case executing the MEX file will
cause MATLAB to crash. But there's only one way to find out: the hard
way.

**Final installation steps.** Installing the mex file is as simple as
placing it wherever MATLAB or Octave can find it; that is, in its
"path". There are many ways to do this. Whatever approach you take, note
that the mex file does not contain the usage documentation, which is
provided in lbfgsb.m. So wherever you place lbfgsb.mex, make sure that
you include lbfgsb.m in the same directory.

Perhaps the simplest approach to copy lbfgsb.mex and lbfgsb.m to a a
specific folder. Alternatively, 

    make install

will create the .mex and .m files in the INSTALLDIR directory, and
this directory will be automatically created if it does not
exist.

Once you are in MATLAB or Octave, add the installation directory to
the MATLAB or Octave path can be done using the command **addpath**.
Note that once you exit MATLAB or Octave the updated path is not
preserved, and you will need to use addpath in every new session.  You
can add the addpath command to ~/.octaverc or ~/matlab/startup. In
this way, the path will automatically be updated every time you start
a new session.

###A brief tutorial

I've written a short script [examplehs038.m](src/examplehs038.m) which
demonstrates how to call **lbfgsb** for solving a very small
optimization problem, the Hock & Schittkowski test problem #38 **(Note
 #3)**. The optimization problem has 4 variables. They are bounded from
below by -10 and from above by +10. We've set the starting point to
(-3, -1, -3, -1). We've written two MATLAB functions for computing the
objective function and the gradient of the objective function at a
given point. These functions can be found in the files
[computeObjectiveHS038.m](src/computeObjectiveHS038.m) and
[computeGradientHS038.m](src/computeGradientHS038.m). If you run the
script, the solver should very quickly progress toward the point (1,
1, 1, 1). The objective is 0 at this point.

The basic MATLAB function call is

    lbfgsb(x0,lb,ub,objfunc,gradfunc)

The first input argument **x0** declares the starting point for the
solver. The second and third input arguments define the lower and
upper bounds on the variables, respectively. If an entry of **lb** is
set to **-Inf**, then this is equivalent to no lower bound for that
entry (the same applies to upper bounds, except that it must be
**+Inf** instead). The next two input arguments must be the names of
MATLAB routines (M-files). The first routine must calculate the
objective function at the current point. The output must be a scalar
representing the objective evaluated at the current point. The second
function must compute the gradient of the objective at the current
point. The input is the same as **objfunc**, but it must return as
many outputs as there are inputs. For a complete guide on using the
MATLAB interface, type **help lbfgsb** in MATLAB.

The script [exampleldaimages.m](src/exampleldaimages.m) is a more
complicated example. It uses the L-BFGS-B solver to compute some
statistics of a posterior distribution for a text analysis model
called latent Dirichlet allication (LDA). The model is used to
estimate the topics for a collection of documents. The optimization
problem is to compute an approximation to the posterior, since it is
not possible to compute the posterior exactly. In this case, the
optimization problem corresponds to minimizing the distance between
the true posterior distribution and an approximate distribution. For
more information, see the paper of Blei, Ng and Jordan, 2003 **(Note
 #4)**. In order to run this example, you will need to install the
[lightspeed](http://research.microsoft.com/~minka/software/lightspeed)
for MATLAB developed by Tom Minka at Microsoft Research.

The script starts by generating a bunch of canonical images that
represent topics; if a pixel is white then a word at that position is
more likely to appear in that topic. Note that the order of the pixels
in the image is unimportant, since LDA is a &quot;bag of words&quot;
model. In another figure we show a small sample of the data generated
from the topics. Each image represents a document, and the pixels in
the image depict the word proportions in that document **(Note #5)**.

Next, the script runs the L-BFGS-B solver to find a local minimum to
the variational objective. Buried in the M-file [mflda.m](src/mflda.m)
is the call to the solver:

    [xi gamma Phi] = lbfgsb({xi gamma Phi},...
        {repmat(lb,W,K)  repmat(lb,K,D)  repmat(lb,K,M)},...
        {repmat(inf,W,K) repmat(inf,K,D) repmat(inf,K,M)},...
        'computeObjectiveMFLDA','computeGradientMFLDA',...
        {nu eta L w lnZconst},'callbackMFLDA',...
        'm',m_lbfgs,'factr',1e-12,'pgtol',...
        tolerance,'maxiter',maxiter);

This example will give us the opportunity to demonstrate some of the
more complicated aspects of the our MATLAB interface. The first input
to **lbfgsb** sets the initial iterate of the optimization
algorithm. Notice that we've passed a cell array with entries that are
matrices. Indeed, it is possible to pass either a matrix or a cell
array. Whatever structure is passed, the other input arguments (and
outputs) must also abide by this structure. The lower and upper bounds
on variables are also cell arrays. All the upper bounds are set to
infinity, which means that the variables are only bounded from
below.

It is instructive to examine the callback routines. They look like:

    f = computeObjectiveMFLDA(xi,gamma,Phi,auxdata)

    [dxi, dgamma, dPhi] = ...
        computeGradientMFLDA(xi,gamma,Phi,auxdata)

Notice that each entry to the cell array is its own input argument to
the callback routines. The same applies for the output arguments of
the gradient callback function.

The sixth input argument is a cell array that contains auxiliary
data. This data will be passed to the MATLAB functions that evaluate
the objective and gradient. The seventh input argument specifies a
callback routine that is called exactly once for every iteration. This
callback routine is most useful for examining the progress of the
solver. The remaining input arguments are label/value pairs that
override some of the default settings of L-BFGS-B.

After converging to a solution, we display topic samples drawn from
the variational approximation. A good solution should come close to
recovering the true topics from which the data was generated (although
they might not be in the same order).

Don't be surprised if it takes many iterations before the optimization
algorithm converges on a local minimum. (In repeated trials, I found
it as many as 2000 iterations to converge to a stationary point of the
objective.) This particular problem demonstrates the limitations of
limited-memory quasi-Newton approximations to the Hessian; the low
storage requirements can come at the cost of slow convergence to the
solution.

###The C++ interface

One of the nice byproducts of writing a MATLAB interface for
L-BFGS-B is that I've ended up with a neat little C++ class that
encapsulates the nuts and bolts of executing the solver. A brief
description of the Program class can be found in the header file
[program.h](src/program.h).

"Program" is an abstract class, which means that its impossible to
instantiate a Program object. This is because the class has member
functions that aren't yet implemented. These are called pure virtual
functions. In order to use the Program class, one needs to define
another class that inherits the Program class and that implements the
pure virtual functions. The class MatlabProgram (in
[matlabprogram.h](src/matlabprogram.h)) is an example of such a class.

The new child class must declare and implement these two
functions:

    virtual double computeObjective (int n, double* x);
    virtual void computeGradient (int n, double* x, double* g);

See the header file for a detailed description of these functions. The
only remaining detail is calling the Program constructor. After that,
it is just a matter of declaring a new object of type MyProgram then
calling the function runSolver.

###Overview of source files

<b>solver.f</b><br>
Fortran 77 source for the L-BFGS-B solver
routines.

<b>program.h<br>
program.cpp<br></b>
Header and source file for the C++ interface.

<b>array.h<br>
arrayofmatrices.h<br>
matlabexception.h<br>
matlabmatrix.h<br>
matlabprogram.h<br>
matlabscalar.h<br>
matlabstring.h<br>
arrayofmatrices.cpp<br>
matlabexception.cpp<br>
matlabmatrix.cpp<br>
matlabprogram.cpp<br>
matlabscalar.cpp<br>
matlabstring.cpp<br>
lbfgsb.cpp<br></b>
Header and source files for the MATLAB interface.

<b>lbfgsb.m</b><br>
Usage instructions for the MEX file.

<b>examplehs038.m</b><br>
MATLAB script for solving the H&S test example #38.

<b>computeObjectiveHS038.m<br>
computeGradientHS038.m<br>
genericcallback.m<br></b>
MATLAB functions used by the script [examplehs038.m](src/examplehs038.m).

<b>exampleldaimages.m</b><br>
MATLAB script that generates synthetic documents and topics, then
computes a mean field variational approximation to the posterior
distribution of the latent Dirichlet allocation model, then displays
the result.

<b>mflda.m</b><br>
This function computes a mean field variational approximation to the
posterior distribution of the latent Dirichlet allocation model by
minimizing the distance between the variational distribution and the
true distribution. It uses the L-BFGS-B to minimize the objective
function subject to the bound constraints. The objective function also
acts as a lower bound on the logarithm of the denominator that appears
from the application of Bayes rule **(Note #4)**.

<b>callbackMFLDA.m<br>
computeObjectiveMFLDA.m<br>
computeGradientMFLDA.m<br>
computelnZ.m<br>
computelnZconst.m<br>
computem.m<br>
dirichletrnd.m<br>
gammarnd.m<br>
createsyntheticdata.m<br>
mfldainit.m<br></b>
Some functions used by [exampleldaimages.m](src/exampleldaimages.m) and
[mflda.m](src/mflda.m).

###Credits

MATLAB interface Developed by<br>
[Peter Carbonetto](http://www.cs.ubc.ca/spider/pcarbo)<br>
Dept. of Human Genetics<br>
University of Chicago

Support for Octave contributed by<br>
[Jose Vallet](http://github.com/josombio)<br>
School of Electrical Engineering<br>
Aalto University

Thank you to Ciyou Zhu, Richard Byrd, Jorge Nocedal and Jose Luis for
making their L-BFGS-B Fortran code available, and to Hannes Nickisch
for ...

###Notes

1. Ciyou Zhu, Richard H. Byrd, Peihuang Lu and Jorge Nocedal
(1997). Algorithm 778: L-BFGS-B: Fortran subroutines for large-scale
bound-constrained optimization. *ACM Transactions on Mathematical
Software* 23(4): 550-560.

2. [Another MATLAB interface](http://www.cs.toronto.edu/~liam/software.shtml)
for the L-BFGS-B routines has been developed by Liam Stewart at the
University of Toronto.

3. Willi Hock and Klaus Schittkowski (1981). *Test Examples for
Nonlinear Programming Codes.* Lecture Notes in Economics and
Mathematical Systems, Vol. 187, Springer-Verlag.

4. David M. Blei, Andrew Y. Ng, Michael I. Jordan (2003). Latent
Dirichlet allocation. *Journal of Machine Learning Research*
Vol3: 993-1022.

5. Tom L. Griffiths and Mark Steyvers (2004). Finding scientific
topics. *Proceedings of the National Academy of Sciences* 101:
5228-5235.
