#include "matlabexception.h"

// Function definitions for class MatlabException.
// -----------------------------------------------------------------
MatlabException::MatlabException (const char* message) throw()
  : exception() { 
  this->message = message;
}

MatlabException::MatlabException (const MatlabException& source) throw() 
  : exception() {
  message = source.message;
}

MatlabException& MatlabException::operator= (const MatlabException& source) 
{ 
  message = source.message; 
  return *this;
}
