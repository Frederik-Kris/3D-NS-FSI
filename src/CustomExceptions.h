/*
 * CustomExceptions.h
 *
 *  Created on: Oct 13, 2022
 *      Author: frederk
 */

#ifndef SRC_CUSTOMEXCEPTIONS_H_
#define SRC_CUSTOMEXCEPTIONS_H_

#include "includes_and_names.h"

template <typename T>
class UnexpectedEnumValueException : public std::exception
{
public:
	UnexpectedEnumValueException( T value )
    {
		std::ostringstream whatStream;
		whatStream << "Value " << value << " of enum " << typeid( T ).name() << " is not supported";
		whatString = whatStream.str();
    }
	const char* what() override
	{
		return whatString.c_str();
	}
private:
	std::string whatString;
};



#endif /* SRC_CUSTOMEXCEPTIONS_H_ */
