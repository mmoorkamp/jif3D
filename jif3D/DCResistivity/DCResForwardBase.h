//============================================================================
// Name        : DCResForwardBase.h
// Author      : Zhanjie Shi and Richard.W Hobbs
// Version     : April 2014
// Copyright   : 2014, Zhanjie Shi and Richard.W Hobbs
//============================================================================

#ifndef DCRESFORWARDBASE_H_
#define DCRESFORWARDBASE_H_

#include <cstdlib>

#include <vector>
#include <iostream>

#include "../Global/Serialization.h"
#include "../Global/VecMat.h"
#include "../Global/Jif3DGlobal.h"

namespace jif3D
  {
    /** \addtogroup Resistivity Parameters classes and functions */
    /***\Referring to RESINVM3D (A Matlab 3D resistivity inversion package, developed by Adam Pidlisecky et.al., 2007, Geophysics, vol72(2), p.H1-H10), we translated Matlab forward code
     * into C++ code. The formulation of forward problem is d(m)=QA'(m)q, A(m)=DS(rho)G. Thereinto, D and G are matrices representing 3D divergence and gradient operators related with
     * model size parameter, respectively; S(rho) is a diagonal matrix related with model size and resistivity values rho; q is a matrix related with the locations of the source electrodes;
     * Q is a projection matrix used to collect data from potential volume, i.e. the matrix related with the locations of receiver electrodes.
     * The implementation method of forward modelling in Matlab program is as follows. Firstly, generate the divergence matrix D, gradient matrix G, diagonal matrix S, source term q and
     * receiver term Q. Secondly, create forward operator A by multiplying D, S, and G which are sparse matrix. Thirdly, calculate the preconditioner by Incomplete LU Decomposition of A.
     * Note that A was converted to permutation using a symmetric reverse Cuthill-McKee permutation (Matlab built-in symcrm.m, Cuthill-McKee, 1969) before ILU. And then, with the
     * preconditioner, acquire potential volume by solving the linear system A(m)u=q using a preconditoned biconjugate, stabilized gradient algorithm(bicgstb.m)(Saad,1996). Finally,
     * collect forward response data from potential volume through the formulation Q*u.
     * The basic function for gradient of objective function is also implemented here because it will need forward operator, Q term and bicgstb function during calculation of gradient.
     * Then this function can be called for LQDerivative() in DCResistivityCalculator.cpp directly. In addition, we only calculate gradient related with data term though the model term
     * was included in Matlab program.
     */
    /* @{ */

    //! Parameters of the grid structure:
    class J3DEXPORT GRID_STRUCT_RES
      {
    private:
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & dx;
          ar & dy;
          ar & dz;
          ar & nx;
          ar & ny;
          ar & nz;
          ar & rho;
        }
    public:
      double dx; /*!< Cell width vector in m in x-direction*/
      double dy; /*!< Cell width vector in m in y-direction*/
      std::vector<double> dz; /*!< Cell width vector in m in z-direction*/
      size_t nx; /*!< The number of cells in x-direction*/
      size_t ny; /*!< The number of cells in y-direction*/
      size_t nz; /*!< The number of cells in z-direction*/
      /*!< Resistivity parameters*/
      std::vector<double> rho; /*!< resistivity model used for the forward model*/
      GRID_STRUCT_RES() :
          dx(0), dy(0), dz(), nx(0), ny(0), nz(0), rho()
        {
        }
      virtual ~GRID_STRUCT_RES()
        {
        }
      };

    //! Geometry of the  the Sources.
    class J3DEXPORT GEOMETRY_RES
      {
    public:
      //! Read and store independent source position only
      std::vector<double> PosSx; /*!< x-coordinates of the locations of positive source electrode in m*/
      std::vector<double> PosSy; /*!< y-coordinates of the locations of positive source electrode in m*/
      std::vector<double> PosSz; /*!< z-coordinates of the locations of positive source electrode in m*/
      std::vector<double> NegSx; /*!< x-coordinates of the locations of negative source electrode in m*/
      std::vector<double> NegSy; /*!< y-coordinates of the locations of negative source electrode in m*/
      std::vector<double> NegSz; /*!< z-coordinates of the locations of negative source electrode in m*/
      size_t nsource; /*!< Number of independent source pairs*/
      //! Read and store all receiver position of each data and corresponding source number
      std::vector<double> rx1; /*!< x-coordinates of the locations of the first receiver electrode for each data in m*/
      std::vector<double> ry1; /*!< y-coordinates of the locations of the first receiver electrode for each data in m*/
      std::vector<double> rz1; /*!< z-coordinates of the locations of the first receiver electrode for each data in m*/
      std::vector<double> rx2; /*!< x-coordinates of the locations of the second receiver electrode for each data in m*/
      std::vector<double> ry2; /*!< y-coordinates of the locations of the second receiver electrode for each data in m*/
      std::vector<double> rz2; /*!< z-coordinates of the locations of the second receiver electrode for each data in m*/
      std::vector<size_t> sno; /*!< Source index corresponding to each receiver electrodes*/
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & PosSx;
          ar & PosSy;
          ar & PosSz;
          ar & NegSx;
          ar & NegSy;
          ar & NegSz;
          ar & rx1;
          ar & ry1;
          ar & rz1;
          ar & rx2;
          ar & ry2;
          ar & rz2;
          ar & nsource;
          ar & sno;
        }
      GEOMETRY_RES() :
          PosSx(), PosSy(), PosSz(), NegSx(), NegSy(), NegSz(), nsource(0), rx1(), ry1(), rz1(), rx2(), ry2(), rz2(), sno()
        {
        }
      virtual ~GEOMETRY_RES()
        {
        }
      };

    //! The basic DC Resistivity forward function
    jif3D::rvec ResForward(const GEOMETRY_RES &geo, const GRID_STRUCT_RES &grid,
        size_t ndata_res);
    //! The basic DC Resistivity gradient calculation function for objective function
    jif3D::rvec ResGradient(const GEOMETRY_RES &geo, const GRID_STRUCT_RES &grid,
        const jif3D::rvec &wdwmisfit);
    //! Interpolation function used to interpolate source/receiver's position to cell centre
    std::vector<double> Linint(const std::vector<double> &x, const std::vector<double> &y,
        const std::vector<double> &z, double xr, double yr, double zr);

  /* @} */

  }

#endif /* DCRESFORWARDBASE_H_ */

