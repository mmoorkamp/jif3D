//============================================================================
// Name        : JointObjective.cpp
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "JointObjective.h"
#include "../Global/convert.h"
#include <boost/format.hpp>
#include <iostream>
#include <iomanip>

using boost::format;
using boost::io::group;

namespace jiba
  {
    static const std::string MisfitFormat = " %15s ";

    JointObjective::JointObjective(bool Verbose) :
      Objectives(), Weights(), IndividualFits(), Distributor(), PrintMisfit(
          Verbose)
      {

      }

    JointObjective::~JointObjective()
      {

      }

    void JointObjective::ImplDataDifference(const jiba::rvec &Model,
        jiba::rvec &Diff)
      {
        size_t totaldata = 0;
        const size_t nobjective = Objectives.size();
        IndividualFits.resize(nobjective);
        //if we set verbose in the constructor
        //write out a description of the objective functions to the screen
        if (PrintMisfit)
          {
            for (size_t i = 0; i < nobjective; ++i)
              {
                std::cout << format(MisfitFormat) % Names.at(i);
              }
            std::cout << std::endl;
          }
        //go through all the objective function objects and calculate the misfit
        //also count how much data points we have in total
        for (size_t i = 0; i < nobjective; ++i)
          {
            //the call to the Distributor object makes sure that each
            //individual objective functions gets the right type of parameters
            //converted from the inversion parameter vector
            IndividualFits.at(i) = Objectives.at(i)->CalcMisfit(Distributor(
                Model, i));
            if (PrintMisfit)
              {
                std::cout << format(MisfitFormat) % IndividualFits.at(i);
              }
            totaldata += Objectives.at(i)->GetDataDifference().size();
          }
        Diff.resize(totaldata);
        size_t currstart = 0;
        //now form a single misfit vector that combines all the individual misfits
        //the total misfit will be the squared sum, so we have to weight
        //each element by the square root of the weight
        for (size_t i = 0; i < nobjective; ++i)
          {
            const size_t ndata = Objectives.at(i)->GetDataDifference().size();
            ublas::vector_range<jiba::rvec>(Diff, ublas::range(currstart,
                currstart + ndata)) = sqrt(Weights.at(i))
                * Objectives.at(i)->GetDataDifference();
            currstart += ndata;
          }
        if (PrintMisfit)
          {
            std::cout << std::endl;
          }
      }

    jiba::rvec JointObjective::ImplGradient(const jiba::rvec &Model,
        const jiba::rvec &Diff)
      {
        //initialize the gradient vector
        jiba::rvec Gradient(Model.size());
        Gradient.clear();
        const size_t nobjective = Objectives.size();
        IndividualGradNorms.resize(nobjective);
        //add up the contribution to the gradient from each method
        //considering the weighting and the transformation
        //that has been applied to the model parameters
        jiba::rvec CurrGrad(Model.size());
        if (PrintMisfit)
          {
            std::cout << "Individual Grad Norms: \n";
          }
        for (size_t i = 0; i < nobjective; ++i)
          {
            //we calculate the "natural" gradient for each method
            // and then pass it to the corresponding parameter transformation
            //class in the model distributor to account for the effect of those
            //transformations.
            CurrGrad = Distributor.TransformGradient(Model,
                Objectives.at(i)->CalcGradient(Distributor(Model, i)), i);
            //We store the norm of each gradient for information during the inversion
            IndividualGradNorms.at(i) = ublas::norm_2(CurrGrad);
            if (PrintMisfit)
              {
                std::cout << format(MisfitFormat) % IndividualGradNorms.at(i);
              }
            //the total gradient is just the weighted sum of the individual gradients
            Gradient += Weights.at(i) * CurrGrad;
          }
        if (PrintMisfit)
          {
            std::cout << "\n" << std::endl;
          }
        return Gradient;
      }
  }
