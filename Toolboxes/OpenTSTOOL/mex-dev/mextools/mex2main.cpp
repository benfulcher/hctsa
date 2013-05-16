#include "mex2main.h"

set<mxArray_ptr, less<mxArray_ptr> > mxArray::array_table;

mxArray::mxArray(const long m, const long n, const mxComplexity flag) : ptr(NULL), M(m), N(n), 
compl_flag(flag), class_ID(mxDOUBLE_CLASS) 
{
	if (compl_flag == mxREAL) {
		ptr = (double*) malloc(M*N*sizeof(double));	
	}
	
}

mxArray::~mxArray() 
{
	free(ptr);
}

// Clear all mxArrays that are still in array_table.
// All mxArrays that are freed by the user by calling mxDestroyArray
// are already erased from the array_table, so we only free left over
// mxArrays
void mxArray::cleanup() {
	set<mxArray_ptr, less<mxArray_ptr> >::iterator i = array_table.begin();
	
	while (i != array_table.end()) {	
		delete (*i++);
	}
}

// maximal number of input and output arrays for the mex-function call
#define MAX_ARGUMENTS 1024

int main(int argc, char** argv)
{	
	long i;
	
	mxArray* plhs[MAX_ARGUMENTS];
	mxArray* prhs[MAX_ARGUMENTS];
	
	// clear arguments
	for (i=0 ; i < MAX_ARGUMENTS; i++) { plhs[i] = prhs[i] = NULL; }
	
	long nlhs, nrhs;
	
	if (argc < 2) {
		fprintf(stderr, "Mex2main-Driver : Two arguments are needed : filename of input_list and desired number of output arrays (1)\n");
		exit(-1);
	}
	if (argc > 2)
		nlhs = atol(argv[2]);
	else
		nlhs = 1;
	
	if ((nlhs < 0) || (nlhs > MAX_ARGUMENTS)) {
		fprintf(stderr, "Bad number of output arrays requested\n");
		exit(-1);
	}
	
	// parse input_list and load all input arrays

	nrhs = parse_input_list(argv[1], prhs);

	// call the mexFunction
	
	mexFunction(nlhs, plhs, nrhs, (const mxArray**) prhs);
	
	// save output arrays
	
	for (i = 0; i < nlhs; i++) {
		char number[64];
		char* outname = (char*) strdup(argv[0]);
		outname = (char*) realloc(outname, strlen(outname) + 129);
		
		strcat(outname, "_out");
		
		sprintf(number, "%ld", i);
		strcat(outname, number);
		
		strcat(outname, ".dat");
		save_array(plhs[i], outname);
		
		free(outname);
	}
	
	// free all mxArrays that have not been freed previously
	mxArray::cleanup();
	
	return 0;
}

void mexErrMsgTxt(const char* text)
{
	fprintf(stderr, "An Error ocurred : ");
	fprintf(stderr, text);
	fprintf(stderr, "\n");
	exit(-1);
}


#define MAX_LINE_SIZE 32768

mxArray* readmatrix(FILE* fp)
{
	long m, n;
	char buffer[MAX_LINE_SIZE];
	char* token;
	char minus[] = ", \n\r";		// tokens that divide numbers
	long rows = 0;					// counter for number of rows
	long columns = 0;				// counter for number of columns
	
	mxArray* x = NULL;
	
	// first make a dry run to see how many rows and columns we have
	while(fgets(buffer, MAX_LINE_SIZE, fp) != NULL) {
		if (buffer[0] != '#') { 		/* comments start with # */
			if (rows <= 100) {			// check first 100 lines of file to see if number of columns are constant	
				long newcolumns = 0;
				token = strtok(buffer, minus);
				while (token != NULL) {
					newcolumns++;
					token = strtok(NULL, minus);
				}
				if ((columns != 0) && (newcolumns != columns)) {
					fprintf(stderr, "Number of columns must be the same for all rows of input data file\n");
					return NULL;
				}
				columns = newcolumns;	
			}
			rows++;
		}
	} 
	
	if (rows * columns == 0) {
		fprintf(stderr, "There seem to be no data in input file\n");
		return NULL;
	}
	
	// now do real read
	
	rewind(fp);
	
	x = mxCreateDoubleMatrix(rows, columns, mxREAL);
	
	//printf("%d %d   %d %d \n", rows, columns, mxGetM(x), mxGetN(x));
	
	m = 0; 
	while(fgets(buffer, MAX_LINE_SIZE, fp) != NULL) {
		n = 0;
		if (buffer[0] != '#') { 		/* Kommentarzeile mit # einleiten */
			token = strtok(buffer, minus);
			if (m >= rows) {
				fprintf(stderr, "Unexpected change of number of rows while reading input data file\n");
				return NULL;
			}
			if (token != NULL) {
				do 			
				{
					if (n >= columns) {
						fprintf(stderr, "Number of columns must be the same for all rows of input data file (row %ld)\n", m+1);
						return NULL;
					}					
					(mxGetPr(x))[m + n * mxGetM(x)] = atof(token);
					n++;
				} while ((token = strtok(NULL, minus)) != NULL);
				m++;
			}
		}
	} 	
	
	return x;
}	

#define INPUT_LIST_MAX_LINESIZE 8192

long parse_input_list(const char* filename, mxArray** prhs)
{
	long arguments_found = 0;
	
	FILE* fid = fopen(filename, "r");

	if (fid != NULL) {
		char buffer[INPUT_LIST_MAX_LINESIZE];
		char buffer2[INPUT_LIST_MAX_LINESIZE];
		
		char* token;
		char minus[] = "[],;' \n\r";
		char* endchar;
		
		while(fgets(buffer, INPUT_LIST_MAX_LINESIZE-2, fid) != NULL) {
			buffer[INPUT_LIST_MAX_LINESIZE-1] = '\0';
			
			if (buffer[0] != '#') { 		/* Kommentarzeile wird mit # einleiten und ignoriert */
				if (buffer[0] == '[') { 		/* Vektor darf in Form [3.234 -324.34 5] gegeben werden */
					long count = 0;
					
#ifdef DEBUG
					printf("%d : %s\n", strlen(buffer), buffer);
#endif 
					strncpy (buffer2, buffer, INPUT_LIST_MAX_LINESIZE-1);

					token = strtok(buffer+1, minus);

					// first count number of elements		
					while (token != NULL) {
#ifdef DEBUG					
						printf("%s ", token);  
#endif
						if (!strcmp(token, "]")) {
							break;
						} else {
							count++; 
						}

						token = strtok(NULL, minus); 
					}

					if (count > 0) {
						if (buffer2[strlen(buffer2)-1] == '\'') {	// see if vector is transposed
							prhs[arguments_found] = mxCreateDoubleMatrix(count, 1, mxREAL);	
						} else {
							prhs[arguments_found] = mxCreateDoubleMatrix(1, count, mxREAL);	
						}
#ifdef DEBUG	
						printf("Found vector argument of length %ld \n", count);
#endif
						// parse this line again and fill values into vector
				
						count = 0;
						token = strtok(buffer2+1, minus);
						while (token != NULL) {
							//printf("%s ", token); 
							double value = strtod(token, &endchar);	// try to see if we got a number value
							if ((value == 0) && (endchar == token)) {
								continue; //break;
							} else {
								(mxGetPr(prhs[arguments_found]))[count++] = value;
#ifdef DEBUG									 
								printf("%lf ", value);
#endif								
							}

							token = strtok(NULL, minus); 	
						}
#ifdef DEBUG						
						printf("\n");
#endif
						arguments_found++;
					} else {
						fprintf(stderr, "Error parsing vector argument %s \n", buffer2);
						exit(-1);
					}
				} else {
					token = strtok(buffer, minus);

					while (token != NULL) {

						double value = strtod(token, &endchar);	// try to see if we got a number value

						if ((value == 0) && (endchar == token)) {	
							// we got a string value instead of a number
							// this string is interpreted as a file name of a data file

							FILE* fp = fopen(token, "r");
							if (fp != NULL) {
								prhs[arguments_found] = readmatrix(fp);
								if (prhs[arguments_found] == NULL) {
									fprintf(stderr, "Error reading data file %s \n", token);
									exit(-1);
								}
#ifdef DEBUG									
								printf("Found matrix argument\n");
#endif
								fclose(fp);	
							} else {
								fprintf(stderr, "Error opening data file %s \n", token);
								exit(-1);
							}
						} else {	// we got a single number
							prhs[arguments_found] = mxCreateDoubleMatrix(1, 1, mxREAL);		
							(mxGetPr(prhs[arguments_found]))[0] = value;
#ifdef DEBUG								
							printf("Found scalar argument\n");
#endif							
						}
						
						arguments_found++;
						token = strtok(NULL, minus);
					}
				}
			}	
		}	  

		fclose(fid);
	} else {
		fprintf(stderr, "Error opening argument list file %s \n", filename);
		exit(-1);
	}
	
	return arguments_found;
}

void save_array(const mxArray* x, const char* filename) {
	if (x != NULL) {	
		FILE* fid = fopen(filename, "w");

		for (long m = 0; m < mxGetM(x); m++) {
			for (long n = 0; n < mxGetN(x); n++) {
				fprintf(fid, "%lf ", (mxGetPr(x))[m + n * mxGetM(x)]);
			}
			fprintf(fid, "\n");
		}

		fclose(fid);
	} else {
		fprintf(stderr, "Warning: One or more output arguments not assigned\n");
	}
}

