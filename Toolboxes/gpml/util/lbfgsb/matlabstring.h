#ifndef INCLUDE_MATLABSTRING
#define INCLUDE_MATLABSTRING

#include "mex.h"
#include <string>

// Function declarations.
// -----------------------------------------------------------------
// Copy a C-style string (i.e. a null-terminated character array).
char* copystring (const char* source);

// Class MatlabString.
// -----------------------------------------------------------------
// This class encapsulates read-only access to a MATLAB character
// array.
class MatlabString {
public:
  
  // The constructor accepts as input a pointer to a Matlab array,
  // which must be a valid string (array of type CHAR).
  explicit MatlabString (const mxArray* ptr);
  
  // The copy constructor makes a full copy of the source string.
  MatlabString (const MatlabString& source);
  
  // The destructor.
  ~MatlabString();
  
  // Conversion operator for null-terminated string.
  operator const char* () const { return s; };
  
protected:
  char* s;  // The null-terminated string.
  
  // The copy assignment operator is not proper, thus remains
  // protected.
  MatlabString& operator= (const MatlabString& source) 
  { return *this; };
};

#endif
