#ifndef INCLUDE_MATLABEXCEPTION
#define INCLUDE_MATLABEXCEPTION

#include <exception>

// Class MatlabException
// -----------------------------------------------------------------
// This class just makes it easier for me to throw exceptions. Its
// functionality really has nothing to do with MATLAB.
class MatlabException : public std::exception {
public:
  MatlabException (const char* message) throw();
  ~MatlabException()                    throw() { };
  
  // The copy constructor makes a shallow copy.
  MatlabException (const MatlabException& source) throw();
  
  // The copy assignment operator makes a shallow copy as well.
  MatlabException& operator= (const MatlabException& source);
  
  // Return the message string.
  virtual const char* what () const throw() { return message; };
  
private:
  const char* message;  // The error message.
};

#endif
