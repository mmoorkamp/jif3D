//============================================================================
// Name        : JointObjective.cpp
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "JointObjective.h"
#include "../Global/convert.h"
#include <boost/format.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>
#include <iomanip>

using boost::format;
using boost::io::group;

namespace jif3D
  {

    /*! We want to stop the inversion when all of the data misfits are below
     * the expected chi-squared value. If the convergence limit for each individual
     * objective function. If it is positive, this objective function is checked
     * for convergence.
     * @param Objective A joint objective containing one or more individual objective functions
     * @return True when converged, fals otherwise.
     */
    bool CheckConvergence(const jif3D::JointObjective &Objective)
      {
        bool terminate = true;
        for (size_t i = 0; i < Objective.GetIndividualFits().size() - 1; ++i)
          {
            if (Objective.GetIndividualFits().at(i)
                > Objective.GetObjective(i).GetNData())
              {
                terminate = false;
              }
            else
              {
                if (Objective.GetObjective(i).ConvergenceLimit() > 0.0)
                  {
                    std::cout << "Reached target misfit." << std::endl;
                    std::cout << "Objective number: " << i << std::endl;
                    std::cout << "Misfit: " << Objective.GetIndividualFits().at(i)
                        << std::endl;
                    std::cout << "Target: " << Objective.GetObjective(i).GetNData()
                        << std::endl;
                  }
              }
          }
        return terminate;
      }

    static const std::string MisfitFormat = " %15s ";

    JointObjective::JointObjective(bool Verbose) :
        Objectives(), Weights(), IndividualFits(), Distributor(), PrintMisfit(Verbose)
      {

      }

    JointObjective::~JointObjective()
      {

      }

    void JointObjective::ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff)
      {
        const size_t nobjective = Objectives.size();
        IndividualFits.resize(nobjective);
        std::vector<double> times(nobjective);
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
        Diff.resize(nobjective);
        for (size_t i = 0; i < nobjective; ++i)
          {

            boost::posix_time::ptime starttime =
                boost::posix_time::microsec_clock::local_time();
            //the call to the Distributor object makes sure that each
            //individual objective functions gets the right type of parameters
            //converted from the inversion parameter vector
            IndividualFits.at(i) = Objectives.at(i)->CalcMisfit(Distributor(Model, i));
            boost::posix_time::ptime endtime =
                boost::posix_time::microsec_clock::local_time();
            times.at(i) = (endtime - starttime).total_seconds();
            if (PrintMisfit)
              {
                std::cout << format(MisfitFormat) % IndividualFits.at(i);
              }
            Diff(i) = sqrt(Weights.at(i) * IndividualFits.at(i));
          }

        if (PrintMisfit)
          {
            std::cout << std::endl;
            std::cout << "Runtimes: \n";
            for (size_t i = 0; i < nobjective; ++i)
              {
                std::cout << format(MisfitFormat) % times.at(i);
              }
            std::cout << "\n" << std::endl;
          }
      }

    jif3D::rvec JointObjective::ImplGradient(const jif3D::rvec &Model,
        const jif3D::rvec &Diff)
      {
        //initialize the gradient vector
        jif3D::rvec Gradient(Model.size());
        Gradient.clear();
        const size_t nobjective = Objectives.size();
        std::vector<double> times(nobjective);
        IndividualGradNorms.resize(nobjective);
        //add up the contribution to the gradient from each method
        //considering the weighting and the transformation
        //that has been applied to the model parameters
        jif3D::rvec CurrGrad(Model.size());
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
            boost::posix_time::ptime starttime =
                boost::posix_time::microsec_clock::local_time();
            CurrGrad = Distributor.TransformGradient(Model,
                Objectives.at(i)->CalcGradient(Distributor(Model, i)), i);
            boost::posix_time::ptime endtime =
                boost::posix_time::microsec_clock::local_time();
            times.at(i) = (endtime - starttime).total_seconds();
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
            std::cout << "\nRuntimes: \n";
            for (size_t i = 0; i < nobjective; ++i)
              {
                std::cout << format(MisfitFormat) % times.at(i);
              }
            std::cout << "\n" << std::endl;
          }
        return Gradient;
      }
  }
