// this common code block is used to create an ATRIA searcher object either from
// a given points set or load it from an matlab mxArray structure that is given
// as the first function argument to mexFunction
//
// These lines of code have to appear very near to the start of the mexFunction routine :
/* 

	int atria_preprocessing_given = 0;		// flag wheter preprocessing is already given on command line

	ATRIA< point_set<euclidian_distance> >* searcher = 0;	
	// try to see if the first parameter given is an atria structure
	// If this is true, the order of input parameters is shifted by one
	
	// try to see if the first parameter given is an atria structure
	// If this is true, the order of input parameters is shifted by one
	if ((nrhs) && (mxIsStruct(prhs[0]))) {
		atria_preprocessing_given = 1;			
		prhs++; 	// these two lines enable us to use the old argument parsing block without changing it
		nrhs--;
	}

*/
/*

#include "create_searcher.cpp" 

/*

// Also : The following line must be appended to the end of function mexFunction :
/*
	delete searcher;
*/	
	
	point_set<euclidian_distance> points(N,dim, p);	
		
	if (atria_preprocessing_given) {
#ifdef MATLAB_MEX_FILE	
		char* metric = 0;
	
		if (mxIsChar(mxGetField(prhs[-1], 0, "optional"))) {	
    		long buflen = (mxGetM(mxGetField(prhs[-1], 0, "optional")) * mxGetN(mxGetField(prhs[-1], 0, "optional"))) + 1;
 			metric = (char*) mxMalloc(buflen);
        	mxGetString(mxGetField(prhs[-1], 0, "optional"), metric, buflen); 
		}
	
		if ((metric == 0) || (!strncmp("euclidian", metric, strlen(metric)))) {
			searcher = new ATRIA< point_set<euclidian_distance> >(points, prhs[-1]);	// this constructor used the data from the preprocessing 
			if (searcher->geterr()) {
				delete searcher;
				searcher = 0;
			}	
		} else
			printf(" ATRIA preprocessing structure was not created using Euclidian metric; doing preprocessing again\n");
	
		mxFree(metric);	
#endif			 			
	} 
	
	if (searcher == 0) {
		searcher = new ATRIA< point_set<euclidian_distance> >(points);	
	}
	
	if (searcher->geterr()) {
		mexErrMsgTxt("Error preparing searcher");
		return;
	}	 
