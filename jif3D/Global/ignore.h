/*
 * ignore.h
 *
 *  Created on: Feb 11, 2021
 *      Author: max
 */

#ifndef GLOBAL_IGNORE_H_
#define GLOBAL_IGNORE_H_


//Template to shut up compiler warnings about unused variables where we are sure the code has to be like this
//taken from   template<class T> void ignore( const T& ) { }
template<class T> void ignore( const T& ) { }

#endif /* GLOBAL_IGNORE_H_ */
