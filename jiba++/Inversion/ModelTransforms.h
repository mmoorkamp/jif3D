/*
 * ModelTransforms.h
 *
 *  Created on: May 15, 2009
 *      Author: mmoorkamp
 */

#ifndef MODELTRANSFORMS_H_
#define MODELTRANSFORMS_H_

#include "ModelDistributor.h"
namespace jiba
  {
    class LogTransform: public jiba::GeneralModelTransform
      {
    public:
      virtual jiba::rvec Transform(const jiba::rvec &FullModel)
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = exp(FullModel(i));
          return Output;
        }

      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative)
        {

          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = exp(FullModel(i)) * Derivative(i);

          return Output;
        }
      LogTransform()
        {
        }
      virtual ~LogTransform()
        {
        }
      };

    class LogDensityTransform: public jiba::GeneralModelTransform
      {
    public:
      virtual jiba::rvec Transform(const jiba::rvec &FullModel)
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = (exp(-FullModel(i)) + 8500.0) / 5000.0;
          return Output;
        }

      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative)
        {
          const double factor = -1.0 / 5000.0;
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output( i) = factor * exp(-FullModel(i)) * Derivative(i);
              // std::cout <<FullModel(i) << " " << Output(i) << " " << Derivative(i) << std::endl;
            }
          return Output;
        }
      LogDensityTransform()
        {
        }
      virtual ~LogDensityTransform()
        {
        }
      };

    class VelTransform: public jiba::GeneralModelTransform
      {
    public:
      virtual jiba::rvec Transform(const jiba::rvec &FullModel)
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = 1. / FullModel(i);
          return Output;
        }
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative)
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = -1. / (FullModel(i) * FullModel(i)) * Derivative(i);
          return Output;
        }
      VelTransform()
        {
        }
      virtual ~VelTransform()
        {
        }
      };

    class VelDensTransform: public jiba::GeneralModelTransform
      {
    public:
      virtual jiba::rvec Transform(const jiba::rvec &FullModel)
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = (FullModel(i) + 8500.0) / 5000.0;
          return Output;
        }
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative)
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = 1.0 / 5000.0 * Derivative(i);
          return Output;
        }
      VelDensTransform()
        {
        }
      virtual ~VelDensTransform()
        {
        }
      };

  }
#endif /* MODELTRANSFORMS_H_ */
