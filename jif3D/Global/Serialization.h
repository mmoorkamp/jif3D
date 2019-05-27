//============================================================================
// Name        : Serialization.h
// Author      : 11 Sep 2015
// Version     : 
// Copyright   : 2015, mm489
//============================================================================

#ifndef GLOBAL_SERIALIZATION_H_
#define GLOBAL_SERIALIZATION_H_

#ifdef HAVEHPX
#include <hpx/config.hpp>
#include <hpx/runtime/serialization/serialize.hpp>
#include <hpx/runtime/serialization/base_object.hpp>
#include <hpx/runtime/serialization/vector.hpp>
#include <hpx/runtime/serialization/array.hpp>
#include <hpx/runtime/serialization/complex.hpp>
#include <hpx/runtime/serialization/shared_ptr.hpp>
#include <hpx/runtime/serialization/serialization_fwd.hpp>
using hpx::serialization::base_object;
using hpx::serialization::access;
using hpx::serialization::make_array
#else
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
using boost::serialization::base_object;
using boost::serialization::access;
using boost::serialization::make_array;
#endif


#endif /* GLOBAL_SERIALIZATION_H_ */
