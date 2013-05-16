// mex2main.h
// Create a standalone application from an mex file by linking the mex-code
// with this wrapper which gives the basic matlab functionality that 
// the mex-code expects to have

#ifndef MEX2MAIN_H
#define MEX2MAIN_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "set.h"

#define mexPrintf printf

typedef enum {
	mxCELL_CLASS = 1,
	mxSTRUCT_CLASS,
	mxOBJECT_CLASS,
	mxCHAR_CLASS,
	mxSPARSE_CLASS,
	mxDOUBLE_CLASS,
	mxSINGLE_CLASS,
	mxINT8_CLASS,
	mxUINT8_CLASS,
	mxINT16_CLASS,
	mxUINT16_CLASS,
	mxINT32_CLASS,
	mxUINT32_CLASS,
	mxINT64_CLASS,		
	mxUINT64_CLASS,		
	mxUNKNOWN_CLASS = -1
} mxClassID;

typedef enum {
    mxDOUBLE_ARRAY      = 2,    /* start here to align with dispatch */
    mxSPARSE_ARRAY,
    mxCHARACTER_ARRAY,
    mxCELL_ARRAY,
    mxSTRUCTURE_ARRAY,
    mxFLOAT_ARRAY,
    mxINT8_ARRAY,
    mxUINT8_ARRAY,
    mxINT16_ARRAY,
    mxUINT16_ARRAY,
    mxINT32_ARRAY,
    mxUINT32_ARRAY,
    mxINT64_ARRAY,      /* place holder - future enhancements */
    mxUINT64_ARRAY,     /* place holder - future enhancements */
    mxOBJECT_ARRAY,
    mxUNKNOWN_ARRAY,
/*
 *  In compiled programs, the 'type' field will be replaced by one
 *  of the fields below.  The translation between these and the fields
 *  used below will be made as follows:
 *
 *  mccBOOL:        type==SCALAR, 2DDOUBLE, or DOUBLE and logical_flag != 0
 *  mccTEXT:        the type field is mxCHAR
 *  mccINT:     the type field is mxINT
 *  mccIX_B:
 *  mccIX_T:
 *  mccIX_I:        An index type--integers containing bools, text, or ints
 *  mccREAL:        the type field is SCALAR, 2DDOUBLE, or DOUBLE, and pi==0
 *  mccCX:      like mccREAL, but pi != 0
 *  mccCELL:        the type field is mxCELL
 *  mccSTRUCT:      the type field is mxSTRUCT
 *  alloc<0:        corresponds to the mxUNASSIGNED type
 */
    mccMASK =       0x20,
    mccBOOL =       0x20,
    mccTEXT,
    mccINT,
    mccREAL,
    mccCX,
    mccCELL,
    mccSTRUCT,  
    mccSPARSEt,
    mccIX_B,
    mccIX_T,
    mccIX_I
            
} mxArrayType;

typedef enum {
    mxREAL,
    mxCOMPLEX
} mxComplexity;

// struct mxArray_tag {
//     mxName            name;
//     mxArrayType       type;
//     int               scope;
//     mxArray          *link;
//     int               number_of_dims;
//     int               nelements_allocated;
//     struct {
//         unsigned int    scalar_flag : 1;
//         unsigned int    logical_flag : 1;
//         unsigned int    empty_flag : 1;
//         unsigned int    global_flag : 1;
//         unsigned int    on_arraylist_flag : 1;
//         unsigned int    zero_imag_flag : 1;
//         unsigned int    static_flag : 1;
//         unsigned int    colon_flag : 1;
//         unsigned int    creation_flag : 1;
//         unsigned int    private_data_flag : 1;
//         unsigned int    unused : 6;
//         unsigned int    kernel_bits : 8;
//         unsigned int    user_bits : 7;
//         unsigned int    string_flag : 1;
//     }   flags;
//   
//     union {
//         struct {
//             int  m;
//             int  n;
//         }   matrix_dims;
//         int  *dim_array;
//     }   size;
//   
//     union {
//         struct {
//             void  *pdata;
//             void  *pimag_data;
//             int   *ir;
//             int   *jc;
//         }   number_array;
//     
//         struct {
//             mxArray  **cells;
//         }   cell_array;
//     
//         struct {
//             mxArray      **fields;
//             mxName        *field_names;
//             char          *object_classname;
//             int            object_tag;  /* if > 0, structure is object */
//             unsigned int   object_chksum;
//             int            number_of_fields;
//         }   structure_array;
// 
//         struct {
//             mxArray  *pglobal;
//         }   global_array;
//     }   data;
// };

// create a minimal matlab compatible matrix class (up to now, only
// real (== double) matrix are supported)

class mxArray;
typedef mxArray* mxArray_ptr;

class mxArray {
		static set<mxArray_ptr, less<mxArray_ptr> > array_table;	// keep a list of all mxArrays

		double* ptr;	// M by N column major matrix (Matlab=Fortran Style), zero based indexing
		long M,N;
		
		mxComplexity compl_flag;
		mxClassID class_ID;
		
	public: 
		mxArray(const long m, const long n, const mxComplexity flag = mxREAL);	
		~mxArray();

		friend long mxGetM(const mxArray* x);
		friend long mxGetN(const mxArray* x);
		friend double* mxGetPr(const mxArray* x);
		friend mxArray *mxCreateDoubleMatrix(int m, int n, mxComplexity flag);
		friend void mxDestroyArray(mxArray *x);
		
		static void cleanup();	// clear all mxArrays in array_table
};

// give a prototype for the mexFunction call

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]);
void mexErrMsgTxt(const char* text);

inline mxArray *mxCreateDoubleMatrix(int m, int n, mxComplexity flag)
{
	mxArray* x = new mxArray(m, n, flag);
	pair<set<mxArray_ptr, less<mxArray_ptr> >::const_iterator, bool > p;
	
	p = mxArray::array_table.insert(x);
	
	if (!p.second) {
		mexErrMsgTxt("Error creating new matrix : pointer is already in use");
		return 0;
	}

#ifdef DEBUG
	fprintf(stderr, "Table count : %d\n", mxArray::array_table.size());
#endif

	return x;
}

// destroy mxArray structure and erase table entry
inline void mxDestroyArray(mxArray *x)
{
	const long count = mxArray::array_table.erase(x);
	
	fprintf(stderr, "%d arrays erased\n",count);
	fprintf(stderr, "Table count : %d\n", mxArray::array_table.size());

	delete x;	
}

inline long mxGetM(const mxArray* x)
{
	return x->M;
}

inline long mxGetN(const mxArray* x)
{
	return x->N;
}

inline double* mxGetPr(const mxArray* x) 
{
	return x->ptr;
}

inline void *mxMalloc(size_t n) {
	return malloc(n);
}
inline void *mxCalloc(size_t n, size_t	size) {
	return calloc(n, size);
}

inline void mxFree(void *ptr) {
	free(ptr);
}
inline void *mxRealloc(void *ptr, size_t size) {
	return realloc(ptr, size);
}

int mxSetDimensions(mxArray *pa, const int *size, int ndims);
void mxSetPr(mxArray *x, double  *pr);
void *mxGetImagData(const mxArray *x);
void mxSetImagData(mxArray *x, void *pi);

inline int mxGetNumberOfDimensions(const mxArray *array_ptr) {
	return 2;
}

inline const int* mxGetDimensions(const mxArray *array_ptr) {
	int* dims = new int[2];
	
	dims[0] = mxGetM(array_ptr);
	dims[1] = mxGetN(array_ptr);
	
	return dims;
}

// Get imaginary data pointer for numeric array
inline double *mxGetPi(const mxArray *x) {
	return 0;
}

// Set imaginary data pointer for numeric array
void mxSetPi(mxArray *x, double  *pi);

// Determine whether the specified array contains numeric (as opposed 
// to cell or struct) data.
inline bool mxIsNumeric(const mxArray *x) {
	return true;
}

// Determine whether the given array is a cell array.
inline bool mxIsCell(const mxArray *x) {
	return false;
}

// Determine whether the given array contains character data. 
inline bool mxIsChar(const mxArray *x) {
	return false;
}

// Determine whether the given array is a sparse (as opposed to full). 
inline bool mxIsSparse(const mxArray *x) {
	return false;
}

// Determine whether the given array is a structure array.
inline bool mxIsStruct(const mxArray *x) {
	return false;
}

// Determine whether the given array contains complex data.
inline bool mxIsComplex(const mxArray *x) {
	return false;
}

// Determine whether the specified array represents its data as 
// double-precision floating-point numbers.
inline bool mxIsDouble(const mxArray *x) {
	return true;
}

// Determine whether the specified array represents its data as 
// single-precision floating-point numbers.
inline bool mxIsSingle(const mxArray *x) {
	return false;
}

// Determine whether the given array's logical flag is on.
inline bool mxIsLogical(const mxArray *x) {
	return false;
}

// Determine whether the specified array represents its data as 
// signed 8-bit integers.
inline bool mxIsInt8(const mxArray *x) {
	return false;
}

// Determine whether the specified array represents its data as 
// unsigned 8-bit integers.
inline bool mxIsUint8(const mxArray *x) {
	return false;
}

// Determine whether the specified array represents its data as 
// signed 16-bit integers.
inline bool mxIsInt16(const mxArray *x) {
	return false;
}

// Determine whether the specified array represents its data as 
// unsigned 16-bit integers.
inline bool mxIsUint16(const mxArray *x) {
	return false;
}

// Determine whether the specified array represents its data as 
// signed 32-bit integers.
inline bool mxIsInt32(const mxArray *x) {
	return false;
}

// Determine whether the specified array represents its data as 
// unsigned 32-bit integers.
inline bool mxIsUint32(const mxArray *x) {
	return false;
}

// prototypes for argument passing fucntions

long parse_input_list(const char* filename, mxArray** prhs);	// substitute a matlab command call by an file based arguement list 
mxArray* readmatrix(FILE* fp);	// arguments that are given from command line are now read in from ASCII file
void save_array(const mxArray* x, const char* filename); 	// output arguments are written as ASCII file

#endif




