/*
 * MultiSectionTransform.h
 *
 *  Created on: 13.02.2013
 *      Author:  mm489
 */

#ifndef MULTISECTIONTRANSFORM_H_
#define MULTISECTIONTRANSFORM_H_

#include <numeric>
#include <vector>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/shared_ptr.hpp>
#include "GeneralModelTransform.h"

namespace jiba
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

    //! For the cross-gradient we sometimes have to extract to sections from the model vector that are not continuous in memory
    /*! For cases were we want to extract several sections from the generalized model parameters
     * that are not continuous in memory, we cannot simply piece a transform
     * from two different SectionTransform objects. For example for the cross-gradient between seismic and MT
     * we currently store the model in the order slowness, density, conductivity. So we have to extract the
     * first n model parameters and the last n model parameters and form a vector of length 2n that
     * is fed to the cross-gradient function. This class generalizes this problem to an arbitrary number
     * of sections. In principle each section can have a different length n1, n2, ... and the resulting
     * transformed vector will have length n1+n2+...
     */
    class MultiSectionTransform: public jiba::GeneralModelTransform
      {
    private:
      //! The length of the model vector containing the full model
      size_t length;
      //! The indices of the first element of the input vector for each transform
      std::vector<size_t> startindices;
      //! The indices of the last element +1 (C loop convention) of the input vector for each transform
      std::vector<size_t> endindices;
      //! The vector of transforms applied to the different sections
      std::vector<boost::shared_ptr<GeneralModelTransform> > Transforms;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<GeneralModelTransform>(*this);
          ar & length;
          ar & startindices;
          ar & endindices;
          ar & Transforms;
        }
      MultiSectionTransform() :
          length(0)
        {

        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual MultiSectionTransform* clone() const
        {
          return new MultiSectionTransform(*this);
        }
      //! Extract several sections from a long vector and apply corresponding parameter transformations
      /*! This function pieces together an output vector from several
       * sections of the input vector. The number and indices and parameter
       * transformations of these sections are determined from the information
       * added in the Constructor and the AddSection member function. One application
       * for such a transform is the CrossGradient objective function. Here we need
       * velocity and conductivity, for example, from the full model vector that contains
       * velocity, density and conductivity information.
       * @param FullModel The full model vector containing all information in generalized form
       * @return A vector pieced together from different sections of the full model and with appropriate transformations applied
       */
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          using boost::numeric::ublas::subrange;
          //we calculate the total length of the transformed vector from
          //the length of each section
          const size_t outlength = std::accumulate(endindices.begin(), endindices.end(),
              0) - std::accumulate(startindices.begin(), startindices.end(), 0);
          jiba::rvec Result(outlength);

          const size_t nsections = startindices.size();
          //the result will be continuous in memory and always
          //start at 0, each section will just follow after the other
          size_t resultstartindex = 0;
          for (size_t i = 0; i < nsections; ++i)
            {
              const size_t resultendindex = resultstartindex + endindices[i]
                  - startindices[i];
              subrange(Result, resultstartindex, resultendindex) =
                  Transforms[i]->GeneralizedToPhysical(
                      subrange(FullModel, startindices[i], endindices[i]));
              resultstartindex = resultendindex;
            }
          return Result;
        }
      //! Build a model vector of a certain length from a shorter vector
      /*! This is the inverse transform to GeneralizedToPhysical. We therefore
       * have to build a potentially longer model vector from pieces that have been
       * extracted and transformed. The length of the original vector is specified in
       * the constructor and the default value for each element is zero. The information
       * from the input vector FullModel is then copied to the corresponding sections
       * in the output vector.
       * @param FullModel The (potentially shorter) vector containing the information pieced together by GeneralizedToPhysical
       * @return The (potentially longer) output vector where the information in the input vector has been copied to the appropriate section, all other values are zero
       */
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          jiba::rvec Result(length);
          Result.clear();
          const size_t nsections = startindices.size();
          size_t currstart = 0;
          for (size_t i = 0; i < nsections; ++i)
            {
              subrange(Result, startindices[i], endindices[i]) =
                  Transforms[i]->PhysicalToGeneralized(
                      subrange(FullModel, currstart,
                          currstart + endindices[i] - startindices[i]));
              currstart += endindices[i] - startindices[i];
            }
          return Result;
        }
      //! Transform the derivative with respect to the sectioned vector to the full vector
      /*! For the CrossGradient objective function, for example, we calculate the
       * derivative with respect to velocity and conductivity. We have to copy this
       * gradient information to the appropriate sections of the full model vector that
       * contains velocity, density and conductivity. The gradient for parts of the
       * model vector that do not correspond to a section here are set to zero.
       * @param FullModel The model vector containing all information
       * @param Derivative The derivative vector with derivative for the different sections
       * @return The drievative for the full model vector
       */
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          using boost::numeric::ublas::subrange;
          jiba::rvec Result(length);
          Result.clear();
          const size_t nsections = startindices.size();
          size_t currstart = 0;
          for (size_t i = 0; i < nsections; ++i)
            {
              subrange(Result, startindices[i], endindices[i]) =
                  Transforms[i]->Derivative(
                      subrange(FullModel, startindices[i], endindices[i]),
                      subrange(Derivative, currstart,
                          currstart + endindices[i] - startindices[i]));
              currstart += endindices[i] - startindices[i];
            }
          return Result;
        }
      //! Add a new transformation and associated range of indices to the object
      /*! We can add a new section by specifying the startindex and endindex
       * in C/C++ loop convention and a corresponding parameter transform.
       * @param startindex The index of the first element of the section
       * @param endindex The index+1 of the last element of the section
       * @param Trans The corresponding parameter transform
       */
      void AddSection(size_t startindex, size_t endindex,
          boost::shared_ptr<GeneralModelTransform> Trans)
        {
          startindices.push_back(startindex);
          endindices.push_back(endindex);
          Transforms.push_back(Trans);
        }
      //! The constructor needs the length of the full model vector as an argument
      /*! As the transformed model vector is potentially shorter than the input
       * vector, we need to specify the length of the original model vector in
       * the constructor.
       * @param l The length of the full model vector.
       */
      MultiSectionTransform(size_t l) :
          length(l)
        {

        }
      //! We have a second constructor that allows us to specify one initial transform
      MultiSectionTransform(size_t l, size_t startindex, size_t endindex,
          boost::shared_ptr<GeneralModelTransform> Trans) :
          length(l)
        {
          AddSection(startindex, endindex, Trans);
        }
      virtual ~MultiSectionTransform()
        {

        }
      };

  /* @} */
  }

#endif /* MULTISECTIONTRANSFORM_H_ */
