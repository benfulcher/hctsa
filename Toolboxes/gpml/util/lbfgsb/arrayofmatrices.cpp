#ifdef __linux__
  #include "mex.h"
#else
  #include "matrix.h"
#endif

#include "arrayofmatrices.h"

// Function definitions for class ArrayOfMatrices.
// -----------------------------------------------------------------
ArrayOfMatrices::ArrayOfMatrices (const mxArray* ptr) 
  : Array<Matrix*>(getnummatlabmatrices(ptr)) {
  setvalue(0);

  if (mxIsCell(ptr))
      
    // Fill in the entries of the array with the entries of the Matlab
    // cell array.
    for (int i = 0; i < n; i++)
      elems[i] = new Matrix(mxGetCell(ptr,i));
  else
      
    // Fill in the single entry of the array with the Matlab matrix.
    elems[0] = new Matrix(ptr);
}

ArrayOfMatrices::ArrayOfMatrices (const mxArray* ptrs[], int numptrs)
  : Array<Matrix*>(numptrs) {
  setvalue(0);

  for (int i = 0; i < numptrs; i++)
    elems[i] = new Matrix(ptrs[i]);
}
  
ArrayOfMatrices::ArrayOfMatrices (mxArray* ptrs[], 
				  const ArrayOfMatrices& source) 
  : Array<Matrix*>(source.length()) {
  setvalue(0);

  // Create the new Matlab matrices and copy the data from the
  // source.
  for (int i = 0; i < n; i++) {
    int h     = source[i]->height();
    int w     = source[i]->width();
    elems[i]  = new Matrix(ptrs[i],h,w);
    *elems[i] = *source[i];
  }
}

ArrayOfMatrices::ArrayOfMatrices (const ArrayOfMatrices& source) 
  : Array<Matrix*>(source) { }

ArrayOfMatrices::ArrayOfMatrices (double* data, 
				  const ArrayOfMatrices& model) 
  : Array<Matrix*>(model.length()) {
  setvalue(0);
  
  // Repeat for each matrix.
  for (int i = 0; i < n; i++) {
    const Matrix* A = model[i];
    elems[i]        = new Matrix(data,A->height(),A->width());
    data           += A->length();
  }
}

ArrayOfMatrices::~ArrayOfMatrices() {
  for (int i = 0; i < n; i++)
    if (elems[i])
      delete elems[i];
}
  
int ArrayOfMatrices::numelems() const {
  int m = 0;  // The return value.

  // Repeat for each matrix.
  for (int i = 0; i < n; i++)
    m += elems[i]->length();

  return m;
}

ArrayOfMatrices& ArrayOfMatrices::operator= (const ArrayOfMatrices& source) {

  // Copy each matrix.
  for (int i = 0; i < n; i++) 
    *elems[i] = *source.elems[i];

  return *this;
}

int ArrayOfMatrices::getnummatlabmatrices (const mxArray* ptr) {
  int n;  // The return value.
    
  if (mxIsCell(ptr))
    n = mxGetNumberOfElements(ptr);
  else if (mxGetNumberOfDimensions(ptr) == 2 && mxIsDouble(ptr))
    n = 1;
  else throw MatlabException("The Matlab array must either be a cell \
array or a double precision matrix");
  return n;    
}
