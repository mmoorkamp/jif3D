//  Copyright (c) 2014 Thomas Heller
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef HPX_SERIALIZATION_VECTOR_HPP
#define HPX_SERIALIZATION_VECTOR_HPP

#ifdef HAVEHPX
#include <hpx/config.hpp>
#include <hpx/runtime/serialization/serialize.hpp>
#include <hpx/traits/is_bitwise_serializable.hpp>
#endif

#include "VecMat.h"

namespace hpx
  { namespace serialization
      {
        // load vector<T>
        template <typename T, typename Allocator>
        void load_impl(input_archive & ar, boost::numeric::ublas::vector<T,Allocator> & vs,
            boost::mpl::false_)
          {
            // normal load ...
            typedef typename boost::numeric::ublas::vector<T,Allocator>::size_type size_type;
            size_type size;
            ar >> size;//-V128
            if(size == 0) return;

            vs.resize(size);
            for(size_type i = 0; i != size; ++i)
              {
                ar >> vs[i];
              }
          }

        template <typename T, typename Allocator>
        void load_impl(input_archive & ar, boost::numeric::ublas::vector<T,Allocator> & v, boost::mpl::true_)
          {
            if(ar.disable_array_optimization())
              {
                load_impl(ar, v, boost::mpl::false_());
              }
            else
              {
                // bitwise load ...
                typedef typename boost::numeric::ublas::vector<T,Allocator>::value_type value_type;
                typedef typename boost::numeric::ublas::vector<T,Allocator>::size_type size_type;
                size_type size;
                ar >> size;//-V128
                if(size == 0) return;

                v.resize(size);
                load_binary(ar, &v[0], v.size() * sizeof(value_type));
              }
          }

        template <typename T, typename Allocator>
        void serialize(input_archive & ar, boost::numeric::ublas::vector<T,Allocator> & v, unsigned)
          {
            v.clear();
            load_impl(
                ar
                , v
                , typename traits::is_bitwise_serializable<
                typename boost::numeric::ublas::vector<T,Allocator>::value_type
                >::type()
            );
          }

        // save vector<T>
        template <typename T, typename Allocator>
        void save_impl(output_archive & ar, const boost::numeric::ublas::vector<T,Allocator> & vs,
            boost::mpl::false_)
          {
            // normal save ...
            typedef typename boost::numeric::ublas::vector<T,Allocator>::value_type value_type;
            for(const value_type & v : vs)
              {
                ar << v;
              }
          }

        template <typename T, typename Allocator>
        void save_impl(output_archive & ar, const boost::numeric::ublas::vector<T,Allocator> & v,
            boost::mpl::true_)
          {
            if(ar.disable_array_optimization())
              {
                save_impl(ar, v, boost::mpl::false_());
              }
            else
              {
                // bitwise save ...
                typedef typename boost::numeric::ublas::vector<T,Allocator>::value_type value_type;
                save_binary(ar, &v[0], v.size() * sizeof(value_type));
              }
          }

        template <typename T, typename Allocator>
        void serialize(output_archive & ar, const boost::numeric::ublas::vector<T,Allocator> & v, unsigned)
          {
            ar << v.size(); //-V128
            if(v.empty()) return;
            save_impl(
                ar
                , v
                , typename traits::is_bitwise_serializable<
                typename boost::numeric::ublas::vector<T,Allocator>::value_type
                >::type()
            );
          }
      }}

#endif
