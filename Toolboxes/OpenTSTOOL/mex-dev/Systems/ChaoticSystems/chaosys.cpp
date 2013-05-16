// mex -I. -O chaosys.cpp
#include <cmath>
#include "mex.h"

#ifdef MATLAB_MEX_FILE
#undef malloc
#undef realloc
#undef free
#define malloc mxMalloc
#define realloc mxRealloc
#define free mxFree
#define printf mexPrintf
#endif

// include header file with the differential equations
#include "dgl.h"

// include declaration and implementation of the ODE solver, based on ode.f from netlib 
#include "adams_pece_decl.h"
#include "adams_pece_impl.h"

double relerr = 1e-12;
double abserr = 1e-8;

template<class System>
mxArray* integrate(System& system, const double stepwidth, const long steps, const double* initials, const long N_intials)
{
	long i,j;
	
	printf("Number of ODEs : %d\n", system.get_number_equations());
	
	if (N_intials < system.get_number_equations()) {
		mexErrMsgTxt("Number of inital values given does not match number of differetial equations"); 
		return 0;
	}
	
	mxArray* out = mxCreateDoubleMatrix(steps, system.get_number_equations(), mxREAL);	/* allocate output matrix */
    double* const y = (double *) mxGetPr(out);   

	Adams<System> solver(system);
	
	for (j=0; j < system.get_number_equations(); j++) {
		system.v[j] = initials[j]; 		// set inital values
	}
	
	long iflag = 1;
	double t = 0;

	for (i=0; i < steps; i++) {
		double tout = t+stepwidth;
		
		solver.ode(&t, &tout, &relerr, &abserr, &iflag);
	
		if (iflag != 2) {
			switch(labs(iflag)) {
				case 3 :
					printf("Increased error tolerances to finish integration. New relerr and abserr : %lf %lf\n",
						relerr, abserr);
					break;
				case 4 :
					printf("error while integrating, too many steps were needed\n");
					mexErrMsgTxt("exiting");
					return out;
					break;
				case 5 :
					mexErrMsgTxt("equations seem to be stiff, exiting");
					return out;
					break;
				case 6:
					mexErrMsgTxt("invalid input parameters, exiting");
					return out;
					break;
				default:
					mexErrMsgTxt("error while integrating, exiting");
					return out;
					break;
			}

			mexPrintf("error in de: iflag = %d ", iflag);
			mexErrMsgTxt("Exiting");
			return out;
        }
		
		for (j=0; j < system.get_number_equations(); j++) {
			y[i + j * steps] = system.v[j];	
		}
		
		//printf("%f %f %f %f\n", t, system.v[0], system.v[1], system.v[2]);
	}
	
	return out;
}

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{       
	/* check input args */
	
	int extra_parameters_given = 0;
	double* params = 0;
	long number_params = 0;

	int mode = 0;
	
	if (nrhs < 3)
	{
		mexErrMsgTxt("chaosys : Number of samples, stepwidth, vector with inital conditions must be given\n, mode, parameter vector and [relerror abserror] are optional");
		return;
	}
	
	/* handle matrix I/O */
	
	const long steps = (long) *(double *)mxGetPr(prhs[0]);
	const double stepwidth = *(double *)mxGetPr(prhs[1]);
	const double* initials = (double *)mxGetPr(prhs[2]);
	
	if (nrhs > 3) {
		mode = (int) *(double *)mxGetPr(prhs[3]);
	}
	
	if (nrhs > 4) {
		extra_parameters_given = 1;
		params = (double *)mxGetPr(prhs[4]);
		number_params = mxGetM(prhs[4])*mxGetN(prhs[4]);
	}
	if (nrhs > 5) {	
		if (mxGetM(prhs[5])*mxGetN(prhs[5]) < 2) {
			mexErrMsgTxt("Argument must be [relerror abserror]"); 
		}
		relerr = ((double *)mxGetPr(prhs[5]))[0];
		abserr = ((double *)mxGetPr(prhs[5]))[1]; 
	}	
	if (steps < 1) {
		mexErrMsgTxt("Number of steps must be positive");
		return;
	}			
	if (stepwidth <= 0) {
		mexErrMsgTxt("Stepsize must be positive");
		return;
	}	
	
	switch(mode) {
		case 1 : {
			printf("Generalized Chua : Double Scroll\n");
			if (extra_parameters_given) {
				if (number_params < 2) {	
					mexErrMsgTxt("wrong number of parameter values"); 
				}
				DoubleScroll dgl(params[0], params[1]);
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2])); 
			} else {
				DoubleScroll dgl;
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2]));	
			}
			break;
		}
		case 2 : {
			printf("Generalized Chua : Five Scroll\n");
			if (extra_parameters_given) {
				if (number_params < 2) {	
					mexErrMsgTxt("wrong number of parameter values"); 
				}
				FiveScroll dgl(params[0], params[1]);
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2])); 
			} else {
				FiveScroll dgl;
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2]));	
			}
			break;
		}		
		case 3 : {
			printf("Duffing\n");		
			if (extra_parameters_given) {
				if (number_params < 3) {	
					mexErrMsgTxt("wrong number of parameter values"); 
				}
				DuffingDGL dgl(params[0], params[1], params[2]);
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2])); 
			} else {
				DuffingDGL dgl;
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2]));	
			}
			break;
		}
		case 4 : {
			printf("Roessler\n");
			if (extra_parameters_given) {
				if (number_params < 3) {	
					mexErrMsgTxt("wrong number of parameter values"); 
				}
				RoesslerDGL dgl(params[0], params[1], params[2]);
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2])); 
			} else {
				RoesslerDGL dgl;
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2]));	
			}
			break;
		}
		case 5 : {
			printf("Toda Oscillator\n");
			if (extra_parameters_given) {
				if (number_params < 3) {	
					mexErrMsgTxt("wrong number of parameter values"); 
				}
				TodaDGL dgl(params[0], params[1], params[2]);
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2])); 
			} /*else {
				TodaDGL dgl;
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2]));	
			}*/
			break;
		}		
		case 6 : {
			printf("Van der Pol Oscillator\n");
			if (extra_parameters_given) {
				if (number_params < 4) {	
					mexErrMsgTxt("wrong number of parameter values"); 
				}
				VanDerPolDGL dgl(params[0], params[1], params[2], params[3]);
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2])); 
			} /*else {
				VanDerPolDGL dgl;
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2]));	
			}*/
			break;
		}		
		case 7 : {
			printf("Pendelum\n");
			if (extra_parameters_given) {
				if (number_params < 4) {	
					mexErrMsgTxt("wrong number of parameter values"); 
				}
				PendelumDGL dgl(params[0], params[1], params[2], params[3]);
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2])); 
			} /*else {
				PendelumDGL dgl;
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2]));	
			}*/
			break;
		}
		case 8 : {
			printf("BaierSahle\n");
			if (extra_parameters_given) {
				if (number_params < 5) {	
					mexErrMsgTxt("wrong number of parameter values"); 
				}
				if (params[0] < 2) {
					mexErrMsgTxt("Number of OD equations must be at least 2"); 
				}
				BaierSahleDGL dgl((long) params[0], params[1], params[2], params[3], params[4]);
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2])); 
			} else {
				BaierSahleDGL dgl;
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2]));	
			}
			break;
		}		
		    
		case 9 : {
			printf("Colpitts\n");
			if (extra_parameters_given) {
				if (number_params < 7) {	
					mexErrMsgTxt("wrong number of parameter values"); 
				}
				ColpittsDGL dgl(params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2])); 
			} else {
			        mexErrMsgTxt("No default parameters available for Colpitts!");
//				ColpittsDGL dgl;
//				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2]));	
			}
			break;
		}
		// case 8 : MackayGlass etc.
								
		default : {
			printf("Lorenz\n");
			if (extra_parameters_given) {
				if (number_params < 3) {	
					mexErrMsgTxt("wrong number of parameter values"); 
				}
				LorenzDGL dgl(params[0], params[1], params[2]);
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2]));
			} else { 
				LorenzDGL dgl;
				plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2]));
			}
			break;	
		}
	}
}

		/*
		case 5 : {
			printf("KuramotoSivashinski\n");
			if (nrhs > 4) {
				const double* params = (double *)mxGetPr(prhs[4]);
				if (mxGetM(prhs[4])*mxGetN(prhs[4]) >= 5) {
					KuramotoSivashinskiDGL dgl(params[0], params[1], params[2], params[3], params[4]); 
					plhs[0] = integrate(dgl, stepwidth, steps, initials, mxGetM(prhs[2])*mxGetN(prhs[2]));	
					break;
				}
			} 
			mexErrMsgTxt("KuramotoSivashinskiDGL : need paramter vector [N A1 A2 A3 A4]");
			break;
		}*/

