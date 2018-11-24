/*
 * MultiSectionCovariance.h
 *
 *  Created on: 24.11.2018
 *      Author:  mm489
 */

#ifndef MULTISECTIONCOVARIANCE_H_
#define MULTISECTIONCOVARIANCE_H_

#include <numeric>
#include <vector>

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/FatalException.h"
#include "GeneralCovariance.h"
#include <boost/shared_ptr.hpp>

namespace jif3D
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
    class J3DEXPORT MultiSectionCovariance: public jif3D::GeneralCovariance
      {
    private:
      //! The length of the model vector containing the full model
      size_t length;
      //! The indices of the first element of the input vector for each transform
      std::vector<size_t> startindices;
      //! The indices of the last element +1 (C loop convention) of the input vector for each transform
      std::vector<size_t> endindices;
      //! The vector of transforms applied to the different sections
      std::vector<boost::shared_ptr<GeneralCovariance >> Transforms;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<GeneralCovariance>(*this);
          ar & length;
          ar & startindices;
          ar & endindices;
          ar & Transforms;
        }
      MultiSectionCovariance() :
          length(0)
        {

        }
    public:

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
      virtual jif3D::rvec ApplyCovar(const jif3D::rvec &vector) override
        {
          using boost::numeric::ublas::subrange;
          //we calculate the total length of the transformed vector from
          //the length of each section
          const size_t outlength = std::accumulate(endindices.begin(), endindices.end(),
              0) - std::accumulate(startindices.begin(), startindices.end(), 0);
          jif3D::rvec Result(outlength);

          const size_t nsections = startindices.size();
          //the result will be continuous in memory and always
          //start at 0, each section will just follow after the other
          size_t resultstartindex = 0;
          for (size_t i = 0; i < nsections; ++i)
            {
              const size_t resultendindex = resultstartindex + endindices[i]
                  - startindices[i];
              subrange(Result, resultstartindex, resultendindex) =
                  Transforms[i]->ApplyCovar(
                      subrange(vector, startindices[i], endindices[i]));
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
      virtual jif3D::rvec ApplyInvCovar(const jif3D::rvec &vector) override
        {
          jif3D::rvec Result(length, 0.0);
          const size_t nsections = startindices.size();
          size_t currstart = 0;
          for (size_t i = 0; i < nsections; ++i)
            {
              subrange(Result, startindices[i], endindices[i]) =
                  Transforms[i]->ApplyInvCovar(
                      subrange(vector, currstart,
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
          boost::shared_ptr<GeneralCovariance> Trans)
        {
          startindices.push_back(startindex);
          endindices.push_back(endindex);
          Transforms.push_back(Trans);
        }
      //! In some cases we need to change the section startindex and endindex after the initial setup
      /*! When initially setting the indices for the transformation in the joint inversion, we only
       * have information on the core grid size. For regularization and cross-gradient this works well.
       * However, for MT and Gravity we might want additional inversion parameters such as background densities
       * or distortion parameters. In order to keep the separation between setting up the coupling
       * and setting up the different methods intact, we need this function to re-adjust the indices
       * when setting the individual functions.
       * @param section The index of the section we want to change
       * @param startindex The index of the first element of the model vector we want to use for this section
       * @param endindex The index of the last element of the model vector we want to use for this section
       */
      void ChangeSectionIndices(size_t section, size_t startindex, size_t endindex)
        {
          if (section >= startindices.size())
            {
              throw jif3D::FatalException(
                  "Trying to change section that does not exist !", __FILE__, __LINE__);
            }
          startindices.at(section) = startindex;
          endindices.at(section) = endindex;
        }

      void SetLength(size_t l)
        {
          length = l;
        }
      //! The constructor needs the length of the full model vector as an argument
      /*! As the transformed model vector is potentially shorter than the input
       * vector, we need to specify the length of the original model vector in
       * the constructor.
       * @param l The length of the full model vector.
       */
      MultiSectionCovariance(size_t l) :
          length(l)
        {

        }
      //! We have a second constructor that allows us to specify one initial transform
      MultiSectionCovariance(size_t l, size_t startindex, size_t endindex,
          boost::shared_ptr<GeneralCovariance> Trans) :
          length(l)
        {
          AddSection(startindex, endindex, Trans);
        }
      virtual ~MultiSectionCovariance()
        {

        }
      };

  /* @} */
  }

#endif /* MULTISECTIONCOVARIANCE_H_ */
