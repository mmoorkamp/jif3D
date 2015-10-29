/*
 * ChainedTransform.h
 *
 *  Created on: 13.02.2013
 *      Author: mm489
 */

#ifndef CHAINEDTRANSFORM_H_
#define CHAINEDTRANSFORM_H_

#include <vector>
#include "../Global/Serialization.h"
#include <boost/shared_ptr.hpp>
#include "GeneralModelTransform.h"

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

//!We can use this transformation class to chain a number of transformations together
    /*! This class takes a number of transformations and consecutively applies them
     * to the generalized parameters. i.e. if we add the transforms f, g and h in that order,
     * we calculate the physical parameters p from the generalized parameters m as
     * \f$ p = h(g(f(m))) \f$ .
     * For the back transformation to physical
     * parameters the order is reversed and we use the chain rule to calculate
     * the derivatives.
     */
    class ChainedTransform: public jif3D::GeneralModelTransform
      {
    private:
      //! We store pointers to each transform in the chain in this vector
      std::vector<boost::shared_ptr<GeneralModelTransform> > Transforms;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<GeneralModelTransform>(*this);
          ar & Transforms;
        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual ChainedTransform* clone() const
        {
          return new ChainedTransform(*this);
        }
      //! Transform the normalized model parameters back to physical parameters
      virtual jif3D::rvec GeneralizedToPhysical(const jif3D::rvec &FullModel) const
        {
          jif3D::rvec Output(FullModel);

          for (size_t j = 0; j < Transforms.size(); ++j)
            Output = Transforms.at(j)->GeneralizedToPhysical(Output);
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jif3D::rvec PhysicalToGeneralized(const jif3D::rvec &FullModel) const
        {
          jif3D::rvec Output(FullModel);
          for (int i = Transforms.size() - 1; i >= 0; --i)
            Output = Transforms.at(i)->PhysicalToGeneralized(Output);
          return Output;
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jif3D::rvec Derivative(const jif3D::rvec &FullModel,
          const jif3D::rvec &Derivative) const
        {
          jif3D::rvec Output(Derivative);
          jif3D::rvec TransModel(FullModel);
          for (size_t j = 0; j < Transforms.size(); ++j)
            {
              Output = Transforms.at(j)->Derivative(TransModel, Output);
              TransModel = Transforms.at(j)->GeneralizedToPhysical(TransModel);
            }
          return Output;
        }
      //! Add a transform to the back of the chain
      void AppendTransform(boost::shared_ptr<GeneralModelTransform> Trans)
        {
          Transforms.push_back(Trans);
        }
      //! Add a transform to the front of the chain
      void PrependTransform(boost::shared_ptr<GeneralModelTransform> Trans)
        {
          Transforms.insert(Transforms.begin(), Trans);
        }
      ChainedTransform()
        {

        }
      virtual ~ChainedTransform()
        {
        }
      };
  /* @} */
  }

#endif /* CHAINEDTRANSFORM_H_ */
