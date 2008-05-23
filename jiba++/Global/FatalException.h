#ifndef FATALEXCEPTION_H_
#define FATALEXCEPTION_H_
#include <stdexcept>
/*! The FatalException class is thrown when there is a problem that cannot be fixed within the program
 * It is derived from runtime_error to make error message handling easier.
 */

class FatalException: public std::runtime_error
{
public:
	FatalException(const std::string whatString) : std::runtime_error(whatString){};
	virtual ~FatalException() throw(){};
};

#endif /*CFATALEXCEPTION_H_*/
