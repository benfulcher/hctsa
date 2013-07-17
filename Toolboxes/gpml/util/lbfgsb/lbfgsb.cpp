#include "mex.h"
#include "matrix.h"
#include "matlabexception.h"
#include "matlabscalar.h"
#include "matlabstring.h"
#include "matlabmatrix.h"
#include "arrayofmatrices.h"
#include "program.h"
#include "matlabprogram.h"
#include <string.h>
#include <exception>

extern void _main();

// Constants.
// -----------------------------------------------------------------
const int minNumInputArgs = 5;

// Function definitions. 
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[]) 
  try {

    // Check to see if we have the correct number of input and output
    // arguments.
    if (nrhs < minNumInputArgs)
      throw MatlabException("Incorrect number of input arguments");

    // Get the starting point for the variables. This is specified in
    // the first input argument. The variables must be either a single
    // matrix or a cell array of matrices.
    int k = 0;  // The index of the current input argument.
    ArrayOfMatrices x0(prhs[k++]);

    // Create the output, which stores the solution obtained from
    // running IPOPT. There should be as many output arguments as cell
    // entries in X.
    if (nlhs != x0.length())
      throw MatlabException("Incorrect number of output arguments");
    ArrayOfMatrices x(plhs,x0);

    // Load the lower and upper bounds on the variables as
    // ArrayOfMatrices objects. They should have the same structure as
    // the ArrayOfMatrices object "x".
    ArrayOfMatrices lb(prhs[k++]);
    ArrayOfMatrices ub(prhs[k++]);

    // Check to make sure the bounds make sense.
    if (lb != x || ub != x)
      throw MatlabException("Input arguments LB and UB must have the same \
structure as X");

    // Get the Matlab callback functions.
    MatlabString objFunc(prhs[k++]);
    MatlabString gradFunc(prhs[k++]);

    // Get the auxiliary data.
    const mxArray* auxData;
    const mxArray* ptr = prhs[k++];
    if (nrhs > 5) {
      if (mxIsEmpty(ptr))
	auxData = 0;
      else
	auxData = ptr;
    }
    else
      auxData = 0;

    // Get the intermediate callback function.
    MatlabString* iterFunc;
    ptr = prhs[k++];
    if (nrhs > 6) {
      if (mxIsEmpty(ptr))
	iterFunc = 0;
      else
	iterFunc = new MatlabString(ptr);
    }
    else
      iterFunc = 0;

    // Set the options for the L-BFGS algorithm to their defaults.
    int    maxiter = defaultmaxiter;
    int    m       = defaultm;
    double factr   = defaultfactr;
    double pgtol   = defaultpgtol;

    // Process the remaining input arguments, which set options for
    // the IPOPT algorithm.
    while (k < nrhs) {

      // Get the option label from the Matlab input argument.
      MatlabString optionLabel(prhs[k++]);

      if (k < nrhs) {

	// Get the option value from the Matlab input argument.
	MatlabScalar optionValue(prhs[k++]);
	double       value = optionValue;

	if (!strcmp(optionLabel,"maxiter"))
	  maxiter = (int) value;
	else if (!strcmp(optionLabel,"m"))
	  m = (int) value;
	else if (!strcmp(optionLabel,"factr"))
	  factr = value / mxGetEps();
	else if (!strcmp(optionLabel,"pgtol"))
	  pgtol = value;
	else {
	  if (iterFunc) delete iterFunc;
	  throw MatlabException("Nonexistent option");
	}
      }
    }    

    // Create a new instance of the optimization problem.
    x = x0;
    MatlabProgram program(x,lb,ub,&objFunc,&gradFunc,iterFunc,
			  (mxArray*) auxData,m,maxiter,factr,pgtol);    

    // Run the L-BFGS-B solver.
    SolverExitStatus exitStatus = program.runSolver();
    if (exitStatus == abnormalTermination) {
      if (iterFunc) delete iterFunc;
      throw MatlabException("Solver unable to satisfy convergence \
criteria due to abnormal termination");
    }
    else if (exitStatus == errorOnInput) {
      if (iterFunc) delete iterFunc;
      throw MatlabException("Invalid inputs to L-BFGS routine");
    }

    // Free the dynamically allocated memory.
    if (iterFunc)
      delete iterFunc;

  } catch (std::exception& error) {
    mexErrMsgTxt(error.what());
  }

