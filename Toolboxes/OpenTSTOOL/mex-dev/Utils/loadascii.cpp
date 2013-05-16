#include "mextools/mextools.h"

#define MAX_LINE_SIZE 65536

mxArray* readmatrix(FILE* fp)
{
     using namespace std;
	long m, n;
	char buffer[MAX_LINE_SIZE];
	char* token;
	char minus[] = ";, \n\r";		// tokens that divide numbers
	long rows = 0;					// counter for number of rows
	
	mxArray* x = NULL;
	
	// first make a dry run 
	while(fgets(buffer, MAX_LINE_SIZE, fp) != NULL) {
		if (buffer[0] != '#') { 		/* comments start with # */	
				token = strtok(buffer, minus);
				while (token != NULL) {
					rows++;
					token = strtok(NULL, minus);
				}
		}
	} 
	
	if (rows == 0) {
		fprintf(stderr, "There seem to be no data in input file\n");
		return NULL;
	}
	
	// now do real read
	
	rewind(fp);
	
	x = mxCreateDoubleMatrix(rows, 1, mxREAL);
		
	m = 0; 
	while(fgets(buffer, MAX_LINE_SIZE, fp) != NULL) {
		if (buffer[0] != '#') { 		/* Kommentarzeile mit # einleiten */
			token = strtok(buffer, minus);
			if (token != NULL) {
				do 			
				{	
					(mxGetPr(x))[m++] = atof(token);
				} while ((token = strtok(NULL, minus)) != NULL);
			}
		}
	} 	
	
	return x;
}	


void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{		
	if (nrhs < 1)
	{
		mexErrMsgTxt("filename must be given");
		return;
	}
	
	long buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0]) * sizeof(mxChar)) + 8;
	char* filename = new char[buflen];

	mxGetString(prhs[0], filename, buflen);

	FILE* fid = fopen(filename, "r");

	if (fid != NULL) {
		plhs[0] = readmatrix(fid);
		fclose(fid);
	} else {
		mexPrintf("Error opening argument list file %s \n", filename);
		mexErrMsgTxt("");
		return;
	}
	
	delete[] filename; 
}
