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
      private:
          	  jiba::rvec Reference;
    public:
      virtual jiba::rvec Transform(const jiba::rvec &FullModel)
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = exp(FullModel(i))*Reference(i);
          return Output;
        }

      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative)
        {

          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = Reference(i) * exp(FullModel(i)) * Derivative(i);

          return Output;
        }
      LogTransform(const jiba::rvec &Ref):
          Reference(Ref)
        {
        }
      virtual ~LogTransform()
        {
        }
      };

    class LogDensityTransform: public jiba::GeneralModelTransform
      {
      private:
          	  jiba::rvec Reference;
    public:
      virtual jiba::rvec Transform(const jiba::rvec &FullModel)
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = (1./(exp(FullModel(i))*Reference(i)) + 8500.0) / 5000.0;
          return Output;
        }

      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative)
        {
          const double factor = -1.0 / 5000.0;
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output( i) = factor/(Reference(i) * exp(FullModel(i))) * Derivative(i);
              std::cout <<FullModel(i) << " " << Reference(i) << " "<< Output(i) << " " << Derivative(i) << std::endl;
            }
          return Output;
        }
      LogDensityTransform(const jiba::rvec &Ref):
          Reference(Ref)
        {
        }
      virtual ~LogDensityTransform()
        {
        }
      };

    class VelTransform: public jiba::GeneralModelTransform
      {
      private:
    	  jiba::rvec Reference;
    public:
      virtual jiba::rvec Transform(const jiba::rvec &FullModel)
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = FullModel(i) * Reference(i);
          return Output;
        }
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative)
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = Reference(i) * Derivative(i);
          return Output;
        }
      VelTransform(const jiba::rvec &Ref):
      Reference(Ref)
        {
        }
      virtual ~VelTransform()
        {
        }
      };

    class VelDensTransform: public jiba::GeneralModelTransform
      {
      private:
          	  jiba::rvec Reference;
    public:
      virtual jiba::rvec Transform(const jiba::rvec &FullModel)
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
          {
            Output( i) = ( 1./(Reference(i)*FullModel(i)) + 8500.0) / 5000.0;
          }
          return Output;
        }
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative)
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
          {
            Output( i) = - 1.0 / (Reference(i)*5000.0
            		* FullModel(i) * FullModel(i))* Derivative(i);
          }
          return Output;
        }
      VelDensTransform(const jiba::rvec &Ref):
    	  Reference(Ref)
        {
        }
      virtual ~VelDensTransform()
        {
        }
      };

  }
#endif /* MODELTRANSFORMS_H_ */
