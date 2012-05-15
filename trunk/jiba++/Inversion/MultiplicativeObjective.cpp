//============================================================================
// Name        : MultiplicativeObjective.cpp
// Author      : 21 Mar 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include "MultiplicativeObjective.h"
#include <boost/format.hpp>
#include <iostream>
#include <iomanip>

using boost::format;
using boost::io::group;

namespace jiba
  {

    static const std::string MisfitFormat = " %15s ";

    MultiplicativeObjective::MultiplicativeObjective(bool Verbose) :
        JointObjective(Verbose)
      {
        // TODO Auto-generated constructor stub

      }

    MultiplicativeObjective::~MultiplicativeObjective()
      {
        // TODO Auto-generated destructor stub
      }

    //the implementation of the misfit calculation
    void MultiplicativeObjective::ImplDataDifference(const jiba::rvec &Model,
        jiba::rvec &Diff)
      {
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
        Diff.resize(1);
        Diff(0) = 1.0;
        //go through all the objective function objects and calculate the misfit
        for (size_t i = 0; i < nobjective; ++i)
          {
            //the call to the Distributor object makes sure that each
            //individual objective functions gets the right type of parameters
            //converted from the inversion parameter vector
            IndividualFits.at(i) = std::max(Objectives.at(i)->CalcMisfit(Distributor(Model, i)),1e-10);
            if (PrintMisfit)
              {
                std::cout << format(MisfitFormat) % IndividualFits.at(i);
              }
            Diff(0) *= IndividualFits.at(i);
          }
        if (PrintMisfit)
          {
            std::cout << std::endl;
          }
        Diff(0) = sqrt(Diff(0));
      }

    //the implementation of the gradient calculation
    jiba::rvec MultiplicativeObjective::ImplGradient(const jiba::rvec &Model,
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
            for (size_t j = 0; j < nobjective; ++j)
              {
                if (j != i)
                  {
                    CurrGrad *= IndividualFits.at(j);
                  }
              }
            Gradient += CurrGrad;
          }
        if (PrintMisfit)
          {
            std::cout << "\n" << std::endl;
          }
        return Gradient;
      }
  } /* namespace jiba */
