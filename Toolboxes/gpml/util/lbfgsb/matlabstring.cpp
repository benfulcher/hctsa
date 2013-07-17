#include "matlabstring.h"
#include "matlabexception.h"

// Function definitions.
// -----------------------------------------------------------------
char* copystring (const char* source) {
  int   n    = strlen(source);  // The length of the string.
  char* dest = new char[n+1];   // The return value.
  strcpy(dest,source);
  return dest;
}

// Function definitions for class MatlabString.
// -----------------------------------------------------------------
MatlabString::MatlabString (const mxArray* ptr) {
  s = 0;

  // Check to make sure the Matlab array is a string.
  if (!mxIsChar(ptr))
    throw MatlabException("Matlab array must be a string (of type CHAR)");
  
  // Get the string passed as a Matlab array.
  s = mxArrayToString(ptr);
  if (s == 0)
    throw MatlabException("Unable to obtain string from Matlab array");
}

MatlabString::MatlabString (const MatlabString& source) {
  s = 0;
  
  // Copy the source string.
  s = copystring(source.s);
}

MatlabString::~MatlabString() { 
  if (s) 
    mxFree(s); 
}

