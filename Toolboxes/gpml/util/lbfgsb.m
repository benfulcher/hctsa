% LBFGSB  Call the nonlinear bound-constrained solver that uses
%         limited-memory BFGS quasi-Newton updates.
%
%   The basic function call is
%   
%     LBFGSB(x0,lb,ub,objfunc,gradfunc)
%
%   The first input argument x0 is either a matrix or a cell array of
%   matrices. It declares the starting point for the solver.
%
%   The second and third input arguments lb and ub must be of the exact same
%   structure as x0. They declare the lower and upper bounds on the
%   variables, respectively. Set an entry to Inf to declare no bound. 
%
%   The next two input arguments must be the names of MATLAB routines
%   (M-files):
%
%     objfunc         Calculates the objective function at the current
%                     point. The routine must accept as many inputs as cell
%                     entry in x0 (or one input if x0 is a matrix). The
%                     output must be a scalar representing the objective
%                     evaluated at the current point.
%         
%     gradfunc        Computes the gradient of the objective at the current
%                     point. The input is the same as objfunc, but it must
%                     return as many outputs as there are inputs, and each
%                     of the outputs have the same matrix structure as its
%                     corresponding input.
%
%   Optionally, one may choose to pass additional auxiliary data to the
%   MATLAB callback routines listed above through the function call
%   LBFGSB(...,auxdata). If auxdata is the empty matrix, no extra
%   information is passed. It is important to observe that the auxiliary
%   data MAY NOT CHANGE through the course of the L-BFGS optimization! The
%   auxiliary data keep the same values as they possessed in the initial
%   call. If you need variables that change over time, you may want to
%   consider global variables (type HELP GLOBAL).
%
%   LBFGSB(...,auxdata,iterfunc) specifies an additional callback routine
%   which is called once per algorithm iteration. The callback routine must
%   take the form ITERFUNC(T,F,X1,X2,...,AUXDATA). T is the current
%   iteration of the algorithm. F is the current value of the objective. The
%   inputs X1, X2,... are exactly the same as the ones passed to objfunc and
%   gradfunc. Finally, extra information may be passed through the input
%   AUXDATA. No outputs are expected from iterfunc. If iterfunc is the empty
%   string, no routine is called.
%
%   These inputs may be followed by parameter/value pairs to modify the
%   default algorithm options. The algorithm option 'maxiter' sets the
%   maximum number of iterations. The rest of the options 'm', 'factr' and
%   'pgtol' are described in the L-BFGS-B documentation. Note that in this
%   implementation we divide factr by the machine precision, so that the
%   convergence condition is satisfied when it is less than factr, not less
%   than factr times the machine precision.

fprintf('You need to compile the LBFGSB code located in util/lbfgsb/.\n')
fprintf('For Linux systems, you can use util/lbfgsb/Makefile.\n')
fprintf('Please, see also the text in doc/README, section 2.\n')