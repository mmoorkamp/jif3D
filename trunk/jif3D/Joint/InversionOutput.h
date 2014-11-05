//============================================================================
// Name        : InversionOutput.h
// Author      : 25 Nov 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#ifndef INVERSIONOUTPUT_H_
#define INVERSIONOUTPUT_H_

#include "../Inversion/JointObjective.h"
#include "../ModelTransforms/GeneralModelTransform.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

namespace jif3D
  {
    /** \addtogroup joint Joint inversion routines */
    /* @{ */
//! Use a parameter transform to translate the model vector to a model object and write the model to files for plotting and saving
    template<class ModelType>
    void SaveModel(const jif3D::rvec &InvModel,
        const jif3D::GeneralModelTransform &Transform, ModelType &ModelObject,
        const std::string &filename)
      {
        jif3D::rvec TransModel = Transform.GeneralizedToPhysical(InvModel);
        assert(TransModel.size() >= ModelObject.GetNModelElements());
        std::copy(TransModel.begin(),
            TransModel.begin() + ModelObject.GetNModelElements(),
            ModelObject.SetData().origin());
        ModelObject.WriteVTK(filename + ".vtk");
        ModelObject.WriteNetCDF(filename + ".nc");
      }

//! Store the current misfit for all individual objectives with appropriate formating in an output stream
    void StoreMisfit(std::ofstream &misfitfile, const size_t iteration,
        const double Misfit, const jif3D::JointObjective &Objective)
      {
        misfitfile << std::setw(5) << iteration << " " << std::setw(15) << Misfit << " ";
        for (size_t i = 0; i < Objective.GetIndividualFits().size(); ++i)
          {
            misfitfile << std::setw(15) << Objective.GetIndividualFits().at(i) << " ";
          }

        misfitfile << " " << Objective.GetNEval();
        misfitfile << std::endl;
      }

//! Store the current RMS for all data objectives with appropriate formating in an output stream
    void StoreRMS(std::ofstream &rmsfile, const size_t iteration,
        const jif3D::JointObjective &Objective)
      {
        std::vector<double> RMS = Objective.GetRMS();
        rmsfile << std::setw(5) << iteration << " " << std::setw(15) << " ";
        for (size_t i = 0; i < RMS.size(); ++i)
          {
            rmsfile << std::setw(15) << RMS.at(i) << " ";
          }

        rmsfile << " " << Objective.GetNEval();
        rmsfile << std::endl;
      }

//! Store the Weights for all objectives with appropriate formating in an output stream
    void StoreWeights(std::ofstream &rmsfile, const size_t iteration,
        const jif3D::JointObjective &Objective)
      {
        std::vector<double> Weights = Objective.GetWeights();
        rmsfile << std::setw(5) << iteration << " " << std::setw(15) << " ";
        for (size_t i = 0; i < Weights.size(); ++i)
          {
            rmsfile << std::setw(15) << Weights.at(i) << " ";
          }

        rmsfile << " " << Objective.GetNEval();
        rmsfile << std::endl;
      }
  /* @} */
  }
#endif /* INVERSIONOUTPUT_H_ */