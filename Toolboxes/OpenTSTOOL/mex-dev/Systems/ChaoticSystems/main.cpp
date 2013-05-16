#include <iostream>
#include <cmath>

// include header file with the differential equations
#include "dgl.h"

// include declaration and implementation of the ODE solver, based on ode.f from netlib 
#include "adams_pece_decl.h"
#include "adams_pece_impl.h"

int main()
{
	long i, steps;
	double stepwidth;
	long iflag;
	double t, tout, relerr, abserr;	
	
	LorenzDGL dgl;		// create Lorenz DEQ system with default parameters
	Adams<LorenzDGL> solver(dgl);
		
// 	HarmonicDGL dgl(0.1);
// 	Adams<HarmonicDGL> solver(dgl);

	relerr = 1e-12;
	abserr = 1e-8;

	steps = 10000;
	stepwidth = 0.01;
	
	t = 0;
	iflag = 1;
	
 	dgl.v[0] = 0.01; 
 	dgl.v[1] = - 0.02;
 	dgl.v[2] = -0.2;
	
	for (i=0;i<steps;i++) {  
    	tout = t+stepwidth;
    	
// 		if ((i % 10) == 0) {	// slowly vary the parameter of the harmonic oscillator
// 			dgl.increase_w(0.01);	
// 		}
		
    	solver.ode(&t, &tout, &relerr, &abserr, &iflag);

    	if (iflag != 2) {
            printf("error in de: iflag = %d \n", iflag);
            return -1;
    	}

    	printf("%lf %lf %lf %lf\n", t, dgl.v[0], dgl.v[1], dgl.v[2]); 
		//printf("%lf %lf %lf\n", t, dgl.v[0], dgl.v[1]); 
	}	
	
	return 0;
}
