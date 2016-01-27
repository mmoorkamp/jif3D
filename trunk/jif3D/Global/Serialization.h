//============================================================================
// Name        : Serialization.h
// Author      : 11 Sep 2015
// Version     : 
// Copyright   : 2015, mm489
//============================================================================

#ifndef GLOBAL_SERIALIZATION_H_
#define GLOBAL_SERIALIZATION_H_

#ifdef HAVEHPX
#include <hpx/runtime/serialization/serialize.hpp>
#include <hpx/runtime/serialization/base_object.hpp>
#include <hpx/runtime/serialization/vector.hpp>
#include <hpx/runtime/serialization/array.hpp>
#include <hpx/runtime/serialization/serialization_fwd.hpp>
using hpx::serialization::base_object;
using hpx::serialization::access;
#else
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
using boost::serialization::base_object;
using boost::serialization::access;
#endif

namespace jif3D
  {
//#ifdef HAVEHPX
//    //template <typename Base, typename Derived> using base_object = hpx::serialization::base_object<Base, Derived>;
//#else
//    //template<typename Base, typename Derived> using base_object = boost::serialization::base_object<Base, Derived>(Derived &d);
//
//    template<typename ... Args>
//    auto base_object(
//        Args&&... args) -> decltype(boost::serialization::base_object(std::forward<Args>(args)...))
//      {
//        return boost::serialization::base_object(std::forward<Args>(args)...);
//      }
//#endif
  }

#endif /* GLOBAL_SERIALIZATION_H_ */
