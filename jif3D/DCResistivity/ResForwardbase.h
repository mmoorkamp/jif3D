//============================================================================
// Name        : ResForwardbase.h
// Author      : Zhanjie Shi and Richard.W Hobbs
// Version     : Feb 2014
// Copyright   : 2014, Zhanjie Shi and Richard.W Hobbs
//============================================================================

#ifndef RESFORWARDBASE_H_
#define RESFORWARDBASE_H_

#include <vector>
#include <cstdlib>
#include <iostream>
#include "../Global/VecMat.h"
#include <boost/serialization/serialization.hpp>

namespace jif3D
{
/** \addtogroup Resistivity Parameters classes and functions */
/***\Referring to RESINVM3D (A Matlab 3D resistivity inversion package, developed by Adam Pidlisecky et.al., 2007, Geophysics, vol72(2), p.H1-H10), we translated forward code of
 * this program into C++ code.The formulation of forward problem is d(m)=QA'(c)q.     A(c)=DS(c)G, D and G are matrices representing 3D divergence and gradient operators, respectively,
 * S(c) is a diagonal matrix containing the conductivity values c; q is a vector containing the locations of the positive and negtive current sources; Q is a projection matrix for
 * selecting data points from volume, i.e. the matrix related to receivers. In Matlab program, firstly, generating the divergence and gradient operator matrices D and G, then creating
 * S matrix and the three resulting sparse matrix are multiplied together to yield the forward operator matrix A. After that, A was converted to permutation matrix using a symmetric
 * reverse Cuthill-McKee permutation (Matlab built-in symcrm.m) to reshape the problem so that the preconditioner for the iterative solver can be calculated faster(Cuthill-McKee, 1969).
 * Moreover, using incomplete LU to generate the preconditioning matrix PL and PU. With this preconditioner, solving the linear system A(c)u=q using a preconditoned biconjugate,
 * stabilized gradient algorithm(bicgstb.m)(Saad,1996). Finally, construct Q matrix for obtaining data from u, compute modelling data using Qu.
 */
    /* @{ */

    //! Parameters of the grid structure:
    class GRID_STRUCT_RES
      {
    private:
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & nx;
          ar & ny;
          ar & nz;
          ar & dx;
          ar & dy;
          ar & dz;
        }
    public:
      std::size_t nx; /*!< Number of grid cell center in x-direction*/
      std::size_t ny; /*!< Number of grid cell center in y-direction*/
      std::size_t nz; /*!< Number of grid cell center in z-direction*/

      float dx; /*!< Size of the grid cells in m in x-direction*/
      float dy; /*!< Size of the grid cells in m in y-direction*/
      float dz; /*!< Size of the grid cells in m in z-direction*/

      /*!< Conductivity parameters*/
      std::vector<float> rho; /*!< resistivity model used for the forward model*/
      GRID_STRUCT_RES() :
          nx(0), ny(0), nz(0), dx(0), dy(0), dz(0), rho()
        {
        }
      virtual ~GRID_STRUCT_RES()
        {
        }
      };

        //! Geometry of the  the Sources and receivers
        class GEOMETRY_RES
          {
        public:
          std::vector<float> sx1; /*!< x-coordinates of the electrod1 locations of source  in m*/
          std::vector<float> sy1; /*!< y-coordinates of the electrod1 locations of source  in m*/
          std::vector<float> sz1; /*!< z-coordinates of the electrod1 locations of source  in m*/
          std::vector<float> sx2; /*!< x-coordinates of the electrod2 locations of source  in m*/
          std::vector<float> sy2; /*!< y-coordinates of the electrod2 locations of source  in m*/
          std::vector<float> sz2; /*!< z-coordinates of the electrod2 locations of source  in m*/
          std::vector<float> rx1; /*!< x-coordinates of the electrod1 locations of receiver  in m*/
          std::vector<float> ry1; /*!< y-coordinates of the electrod1 locations of receiver  in m*/
          std::vector<float> rz1; /*!< z-coordinates of the electrod1 locations of receiver  in m*/
          std::vector<float> rx2; /*!< x-coordinates of the electrod2 locations of receiver  in m*/
          std::vector<float> ry2; /*!< y-coordinates of the electrod2 locations of receiver  in m*/
          std::vector<float> rz2; /*!< z-coordinates of the electrod2 locations of receiver  in m*/
          std::size_t nsource; /*!< Number of source*/
          std::size_t nreceiver; /*!< Number of receiver for each shot*/
          friend class boost::serialization::access;
          //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
          template<class Archive>
          void serialize(Archive & ar, const unsigned int version)
            {
              ar & sx1;
              ar & sy1;
              ar & sz1;
              ar & sx2;
              ar & sy2;
              ar & sz2;
              ar & rx1;
              ar & ry1;
              ar & rz1;
              ar & rx2;
              ar & ry2;
              ar & rz2;
              ar & nsource;
              ar & nreceiver;
            }
          GEOMETRY_RES() :
              sx1(), sy1(), sz1(), sx2(), sy2(), sz2(), rx1(), ry1(), rz1(),rx2(), ry2(), rz2(), nsource(0), nreceiver(0)
            {
            }
          virtual ~GEOMETRY_RES()
            {
            }
          };

        //! Structure including informations about the electric potential data
        class DATA_STRUCT_RES
          {
        public:
          /*!< Electric potential parameters*/
          std::vector<float> dcal; /*!< Forward data of model */

          friend class boost::serialization::access;
          //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
          template<class Archive>
          void serialize(Archive & ar, const unsigned int version)
            {

              ar & dcal;

            }
          DATA_STRUCT_RES() :
              dcal()
            {
            }
          virtual ~DATA_STRUCT_RES()
            {
            }
          };
        //! The basic Resistivity forward function, uses the Preconditioned biconjugate, stabilized gradient algorithm calculate potential for different source pairs and
        //       then construct measured results for different receiver combines using a projection matrix Q
        int ResForward(const GEOMETRY_RES &geo, const GRID_STRUCT_RES &grid, DATA_STRUCT_RES *data);
        std::vector<float> Interpolate_N(std::vector<float> x, std::vector<float> y, std::vector<float> z, float xr, float yr, float zr);

      /* @} */
      }
    #endif














