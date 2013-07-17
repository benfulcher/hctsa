#ifndef INCLUDE_PROGRAM
#define INCLUDE_PROGRAM

// Type definitions.
// -----------------------------------------------------------------
// This defines the possible results of running the L-BFGS-B solver.
enum SolverExitStatus { 
  success,              // The algorithm has converged to a stationary
			// point or has reached the maximum number of
			// iterations.
  abnormalTermination,  // The algorithm has terminated abnormally and
			// was unable to satisfy the convergence
			// criteria.
  errorOnInput          // The routine has detected an error in the
			// input parameters.
};

// Constants.
// -----------------------------------------------------------------
const int    defaultm          = 5;
const int    defaultmaxiter    = 100;
const double defaultfactr      = 1e7;
const double defaultpgtol      = 1e-5;
const int    defaultprintlevel = -1;

// Function declarations.
// -----------------------------------------------------------------
// This is the L-BFGS-B routine implemented in Fortran 77.
extern "C" void setulb_ (int* n, int* m, double x[], double l[], 
			 double u[], int nbd[], double* f, double g[], 
			 double* factr, double* pgtol, double wa[], 
			 int iwa[], char task[], int* iprint, 
			 char csave[], bool lsave[], int isave[], 
			 double dsave[]);

// Class Program.
// -----------------------------------------------------------------
// This class encapsulates execution of the L-BFGS-B routine for
// solving a nonlinear optimization problem with bound constraints
// using limited-memory approximations to the Hessian.
//
// This is an abstract class since some of the class methods have not
// been implemented (they are "pure virtual" functions). In order to
// use this class, one needs to define a child class and provide
// definitions for the pure virtual methods.
class Program {
public:

  // The first input argument "n" is the number of variables, and "x"
  // must be set to the suggested starting point for the optimization
  // algorithm. See the L-BFGS-B documentation for more information on
  // the inputs to the constructor.
  Program (int n, double* x, double* lb, double* ub, int* btype,
	   int m = defaultm, int maxiter = defaultmaxiter,
	   double factr = defaultfactr, double pgtol = defaultpgtol);

  // This is the same constructor as above, except that the
  // appropriate amount of memory is dynamically allocated for the
  // variables, the bounds, and the bound types. The destructor also
  // makes sure that the memory is properly deallocated, but not so
  // for the other constructor.
  Program (int n, int m = defaultm, int maxiter = defaultmaxiter,
	   double factr = defaultfactr, double pgtol = defaultpgtol);

  // The destructor.
  virtual ~Program();

  // An implementation of this method should return the value of the
  // objective at the current value of the variables "x", where "n" is
  // the number of variables.
  virtual double computeObjective (int n, double* x) = 0;

  // An implementation of this method should fill in the values of the
  // gradient "g" and the current point "x".
  virtual void computeGradient (int n, double* x, double* g) = 0;

  // The child class may optionally override this method which is
  // called once per iteration of the L-BFGS-B routine. The input
  // arguments are the current iteration "t", the current value of the
  // variables "x" and the current value of the objective "f".
  virtual void iterCallback (int t, double* x, double f) { };

  // Run the solver. Upon completion, the solution is stored in "x".
  SolverExitStatus runSolver();

protected:

  // The copy constructor and copy assignment operator are kept
  // private so that they are not used.
  Program            (const Program& source) { };
  Program& operator= (const Program& source) { return *this; };

  int     n;       // The number of variables.
  double* x;       // The current point.
  double* lb;      // The lower bounds.
  double* ub;      // The upper bounds.
  int*    btype;   // The bound types.
private:

  // Execute a single step the L-BFGS-B solver routine.
  void callLBFGS (const char* cmd = 0);

  // Initialize some of the structures used by the L-BFGS-B routine.
  void initStructures();

  bool    owner;   // If true, then the blocks of memory pointed to by
		   // x, lb, ub and btype must be deallocated in the
		   // destructor.
  double  f;       // The value of the objective.
  double* g;       // The value of the gradient.
  int     iprint;  // The print level.
  int     maxiter; // The maximum number of iterations.
  double  pgtol;   // Convergence parameter passed to L-BFGS-B.
  double  factr;   // Convergence parameter passed to L-BFGS-B.
  int     m;       // The number of variable corrections to the 
                   // limited-memory approximation to the Hessian.

  // These are structures used by the L-BFGS-B routine.
  double* wa;
  int*    iwa;
  char    task[60];
  char    csave[60];
  bool    lsave[4];
  int     isave[44];
  double  dsave[29];
};

#endif
