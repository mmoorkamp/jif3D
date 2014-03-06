//============================================================================
// Name        : ResForwardbase.cpp
// Author      : Zhanjie Shi and Richard.W Hobbs
// Version     : Feb 2014
// Copyright   : 2014, Zhanjie Shi and Richard.W Hobbs
//============================================================================

#include "ResForwardbase.h"
#include "../Global/VecMat.h"
#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <Eigen/Eigen/Dense>
#include <Eigen/Eigen/Sparse>
#include <Eigen/Eigen/LU>
#include <Eigen/Eigen/Core>



using namespace Eigen;


namespace jif3D
{
  std::vector<float> Interpolate_N(std::vector<float> x, std::vector<float> y, std::vector<float> z, float xr, float yr, float zr)
  {

	  float minx;
	  size_t imx=0;
	  std::vector<size_t> ind_x(2,0), dx(2,0);
	  minx=fabs(x[0]-xr);
	  for (size_t i=1; i<x.size(); i++)
	  {
		  if (fabs(x[i]-xr)<minx)
			  imx=i;
	  }
	  if ((xr-x[imx]) >= 0.0)
		{
		  ind_x[0]=imx;
		  ind_x[1]=imx+1;
		}
	  else if ((xr-x[imx]) < 0.0)
	  {
		  ind_x[0]=imx-1;
		  ind_x[1]=imx;
	  }
	  dx[0]=xr-x[ind_x[0]];
	  dx[1]=x[ind_x[1]]-xr;


	  float miny;
	  size_t imy=0;
	  std::vector<size_t> ind_y(2,0), dy(2,0);
	  miny=fabs(y[0]-yr);
	  for (size_t i=1; i<y.size(); i++)
	  {
		  if (fabs(y[i]-yr)<miny)
			  imy=i;
	  }
	  if ((yr-y[imy]) >= 0.0)
		{
		  ind_y[0]=imy;
		  ind_y[1]=imy+1;
		}
	  else if ((yr-y[imy]) < 0.0)
	  {
		  ind_y[0]=imy-1;
		  ind_y[1]=imy;
	  }
	  dy[0]=yr-y[ind_y[0]];
	  dy[1]=y[ind_y[1]]-yr;

	  float minz;
	  size_t imz=0;
	  std::vector<size_t> ind_z(2,0), dz(2,0);
	  minz=fabs(z[0]-zr);
	  for (size_t i=1; i<z.size(); i++)
	  {
		  if (fabs(z[i]-zr)<minz)
			  imz=i;
	  }
	  if ((zr-z[imz]) >= 0.0)
		{
		  ind_z[0]=imz;
		  ind_z[1]=imz+1;
		}
	  else if ((zr-z[imz]) < 0.0)
	  {
		  ind_z[0]=imz-1;
		  ind_z[1]=imz;
	  }
	  dz[0]=zr-z[ind_z[0]];
	  dz[1]=z[ind_z[1]]-zr;


	  float Dx=x[ind_x[1]]-x[ind_x[0]];
	  float Dy=y[ind_y[1]]-y[ind_y[0]];
	  float Dz=z[ind_z[1]]-z[ind_z[0]];

	  std::vector<float> v(x.size()*y.size()*z.size(), 0);

      v[ ind_x[0]+ind_y[0]*x.size()+ind_z[0]*x.size()*y.size() ] = (1-dx[0]/Dx)*(1-dy[0]/Dy)*(1-dz[0]/Dz);
      v[ ind_x[0]+ind_y[1]*x.size()+ind_z[0]*x.size()*y.size() ] = (1-dx[0]/Dx)*(1-dy[1]/Dy)*(1-dz[0]/Dz);
      v[ ind_x[1]+ind_y[0]*x.size()+ind_z[0]*x.size()*y.size() ] = (1-dx[1]/Dx)*(1-dy[0]/Dy)*(1-dz[0]/Dz);
      v[ ind_x[1]+ind_y[1]*x.size()+ind_z[0]*x.size()*y.size() ] = (1-dx[1]/Dx)*(1-dy[1]/Dy)*(1-dz[0]/Dz);
      v[ ind_x[0]+ind_y[0]*x.size()+ind_z[1]*x.size()*y.size() ] = (1-dx[0]/Dx)*(1-dy[0]/Dy)*(1-dz[1]/Dz);
      v[ ind_x[0]+ind_y[1]*x.size()+ind_z[1]*x.size()*y.size() ] = (1-dx[0]/Dx)*(1-dy[1]/Dy)*(1-dz[1]/Dz);
      v[ ind_x[1]+ind_y[0]*x.size()+ind_z[1]*x.size()*y.size() ] = (1-dx[1]/Dx)*(1-dy[0]/Dy)*(1-dz[1]/Dz);
      v[ ind_x[1]+ind_y[1]*x.size()+ind_z[1]*x.size()*y.size() ] = (1-dx[1]/Dx)*(1-dy[1]/Dy)*(1-dz[1]/Dz);

      std::vector<float> Q(v);

      return(Q);
  }

  int ResForward(const GEOMETRY_RES &geo, const GRID_STRUCT_RES &grid, DATA_STRUCT_RES *data)
  {
	    std::vector<float> cellwidth_x(grid.nx, grid.dx);
	    std::vector<float> cellwidth_y(grid.ny, grid.dy);
	    std::vector<float> cellwidth_z(grid.nz, grid.dz);

	//*******Generate index and value of non-zero elements of the divergence matrix D which is a sparse matrix.
	std::vector<int> lxd,lyd,lzd,jxd,jyd,jzd; //the index of non-zero elements in sparse matrix D.
    std::vector<float> kxd,kyd,kzd;           //the value of non-zero elements in sparse matrix D.
	 /*****Divergence matrix D[np][nax+nay+naz], np=grid.nx*grid.ny*grid.nz, nax=(grid.nx-1)*grid.ny*grid.nz,
	 * nay=grid.nx*(grid.ny-1)*grid.nz, naz=grid.nx*grid.ny*(grid.nz-1).
	 *for (size_t i=0; i<2*nax; i++), D[lxd[i]][jxd[i]]=kxd[i]
	 *for (size_t i=0; i<2*nay; i++), D[lyd[i]][jyd[i]+nax]=kyd[i]
	 *for (size_t i=0; i<2*naz; i++), D[lzd[i]][jzd[i]+nax+nay]=kzd[i]
	 *****/

	//***generate d/dx
	//*Entries(l,j,k)
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=1; k<grid.nx-1; k++)
    		{
    			lxd.push_back(k+j*grid.nx+i*grid.nx*grid.ny);
    		}
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx-2; k++)
    		{
    			jxd.push_back(k+j*(grid.nx-1)+i*(grid.nx-1)*grid.ny);
    		}
	for (size_t i=0; i<grid.nz; i++)
		for (size_t j=0; j<grid.ny; j++)
			for (size_t k=1; k<grid.nx-1; k++)
			{
				kxd.push_back(-1.0/(cellwidth_x[k]*1.0*1.0));
			}
	//*Entries(l+1,j,k)
	for (size_t i=0; i<(grid.nx-2)*grid.ny*grid.nz; i++)
	{
		lxd.push_back(lxd[i]);
	}
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=1; k<grid.nx-1; k++)
    		{
    			jxd.push_back(k+j*(grid.nx-1)+i*(grid.nx-1)*grid.ny);
    		}
	for (size_t i=0; i<(grid.nx-2)*grid.ny*grid.nz; i++)
	{
		kxd.push_back(-kxd[i]);
	}
	//*BC at x=0
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    	{
    		lxd.push_back(j*grid.nx+i*grid.nx*grid.ny);
    	}
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    	{
    		jxd.push_back(j*(grid.nx-1)+i*(grid.nx-1)*grid.ny);
    	}
	for (size_t i=0; i<grid.nz; i++)
		for (size_t j=0; j<grid.ny; j++)
		{
			kxd.push_back(1.0/(cellwidth_x[0]*1.0*1.0));
		}
	//*BC at x=end
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    	{
    		lxd.push_back(grid.nx-1+j*grid.nx+i*grid.nx*grid.ny);
    	}
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    	{
    		jxd.push_back(grid.nx-2+j*(grid.nx-1)+i*(grid.nx-1)*grid.ny);
    	}
	for (size_t i=0; i<grid.nz; i++)
		for (size_t j=0; j<grid.ny; j++)
		{
			kxd.push_back(-1.0/(cellwidth_x[grid.nx-1]*1.0*1.0));
		}
	//***generate d/dy
	//*Entries(l,j,k)
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=1; j<grid.ny-1; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			lyd.push_back(k+j*grid.nx+i*grid.nx*grid.ny);
    		}
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny-2; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			jyd.push_back(k+j*grid.nx+i*grid.nx*(grid.ny-1));
    		}
	for (size_t i=0; i<grid.nz; i++)
		for (size_t j=1; j<grid.ny-1; j++)
			for (size_t k=0; k<grid.nx; k++)
			{
				kyd.push_back(-1.0/(1.0*cellwidth_y[j]*1.0));
			}
	//*Entries(l+1,j,k)
	for (size_t i=0; i<grid.nx*(grid.ny-2)*grid.nz; i++)
	{
		lyd.push_back(lyd[i]);
	}
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=1; j<grid.ny-1; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			jyd.push_back(k+j*grid.nx+i*grid.nx*(grid.ny-1));
    		}
	for (size_t i=0; i<grid.nx*(grid.ny-2)*grid.nz; i++)
	{
		kyd.push_back(-kyd[i]);
	}
	//*BC on y=0
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t k=0; k<grid.nx; k++)
    	{
    		lyd.push_back(k+i*grid.nx*grid.ny);
    	}
	for (size_t i=0; i<grid.nz; i++)
		for (size_t k=0; k<grid.nx; k++)
		{
			jyd.push_back(k+i*grid.nx*(grid.ny-1));
		}
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t k=0; k<grid.nx; k++)
    	{
    		kyd.push_back(1.0/(1.0*cellwidth_y[0]*1.0));
    	}
	//*BC on y=end
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t k=0; k<grid.nx; k++)
    	{
    		lyd.push_back(k+grid.nx*(grid.ny-1)+i*grid.nx*grid.ny);
    	}
	for (size_t i=0; i<grid.nz; i++)
		for (size_t k=0; k<grid.nx; k++)
		{
			jyd.push_back(k+grid.nx*(grid.ny-2)+i*grid.nx*(grid.ny-1));
		}
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t k=0; k<grid.nx; k++)
    	{
    		kyd.push_back(-1.0/(1.0*cellwidth_y[grid.ny-1]*1.0));
    	}
	//***generate d/dz
	//*Entries(l,j,k)
    for (size_t i=1; i<grid.nz-1; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			lzd.push_back(k+j*grid.nx+i*grid.nx*grid.ny);
    		}
    for (size_t i=0; i<grid.nz-2; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			jzd.push_back(k+j*grid.nx+i*grid.nx*grid.ny);
    		}
	for (size_t i=1; i<grid.nz-1; i++)
		for (size_t j=0; j<grid.ny; j++)
			for (size_t k=0; k<grid.nx; k++)
			{
				kzd.push_back(-1.0/(1.0*1.0*cellwidth_z[i]));
			}
	//*Entries(l+1,j,k)
	for (size_t i=0; i<grid.nx*grid.ny*(grid.nz-2); i++)
	{
		lzd.push_back(lzd[i]);
	}
    for (size_t i=1; i<grid.nz-1; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			jzd.push_back(k+j*grid.nx+i*grid.nx*grid.ny);
    		}
	for (size_t i=0; i<grid.nx*grid.ny*(grid.nz-2); i++)
	{
		kzd.push_back(-kzd[i]);
	}
	//*BC on z=0
	for (size_t j=0; j<grid.ny; j++)
		for (size_t k=0; k<grid.nx; k++)
		{
			lzd.push_back(k+j*grid.nx);
		}
    for (size_t j=0; j<grid.ny; j++)
    	for (size_t k=0; k<grid.nx; k++)
    	{
    		jzd.push_back(k+j*grid.nx);
    	}
	for (size_t j=0; j<grid.ny; j++)
		for (size_t k=0; k<grid.nx; k++)
		{
			kzd.push_back(1.0/(1.0*1.0*cellwidth_z[0]));
		}
	//*BC on z=end
	for (size_t j=0; j<grid.ny; j++)
		for (size_t k=0; k<grid.nx; k++)
		{
			lzd.push_back(k+j*grid.nx+grid.nx*grid.ny*(grid.nz-1));
		}
    for (size_t j=0; j<grid.ny; j++)
    	for (size_t k=0; k<grid.nx; k++)
    	{
    		jzd.push_back(k+j*grid.nx+grid.nx*grid.ny*(grid.nz-2));
    	}
	for (size_t j=0; j<grid.ny; j++)
		for (size_t k=0; k<grid.nx; k++)
		{
			kzd.push_back(-1.0/(1.0*1.0*cellwidth_z[grid.nz-1]));
		}

	//*******Generate index and value of non-zero elements of the gradient matrix G which is a sparse matrix.
	std::vector<int> lxg,lyg,lzg,jxg,jyg,jzg; //the index of non-zero elements in sparse matrix G.
    std::vector<float> kxg,kyg,kzg;           //the value of non-zero elements in sparse matrix G.

	 /*****Gradient matrix G[nax+nay+naz][np], nax=(grid.nx-1)*grid.ny*grid.nz,
	 * nay=grid.nx*(grid.ny-1)*grid.nz, naz=grid.nx*grid.ny*(grid.nz-1),np=grid.nx*grid.ny*grid.nz.
	 *for (size_t i=0; i<2*nax; i++), G[lxg[i]][jxg[i]]=kxg[i];
	 *for (size_t i=0; i<2*nay; i++), G[lyg[i]+nax][jyg[i]]=kyg[i]
	 *for (size_t i=0; i<2*naz; i++), G[lzg[i]+nax+nay][jzg[i]]=kzg[i]
	 *****/

    //***generate d/dx
    //*Entries (l,j,k)
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx-1; k++)
    		{
    			lxg.push_back(k+j*(grid.nx-1)+i*(grid.nx-1)*grid.ny);
    		}
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx-1; k++)
    		{
    			jxg.push_back(k+j*grid.nx+i*grid.nx*grid.ny);
    		}
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx-1; k++)
    		{
    			kxg.push_back(-2.0/(cellwidth_x[k]*1.0*1.0+cellwidth_x[k+1]*1.0*1.0));
    		}
    //*Entries (l+1,j,k)
    for (size_t i=0; i<(grid.nx-1)*grid.ny*grid.nz; i++)
    {
    	lxg.push_back(lxg[i]);
    }
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=1; k<grid.nx; k++)
    		{
    			jxg.push_back(k+j*grid.nx+i*grid.nx*grid.ny);
    		}
    for (size_t i=0; i<(grid.nx-1)*grid.ny*grid.nz; i++)
    {
    	kxg.push_back(-kxg[i]);
    }
    //***generate d/dy
    //*Entries (l,j,k)
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny-1; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			lyg.push_back(k+j*grid.nx+i*grid.nx*(grid.ny-1));
    		}
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny-1; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			jyg.push_back(k+j*grid.nx+i*grid.nx*grid.ny);
    		}
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny-1; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			kyg.push_back(-2.0/(1.0*cellwidth_y[j]*1.0+1.0*cellwidth_y[j+1]*1.0));
    		}
    //*Entries (l+1,j,k)
    for (size_t i=0; i<grid.nx*(grid.ny-1)*grid.nz; i++)
    {
    	lyg.push_back(lyg[i]);
    }
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=1; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			jyg.push_back(k+j*grid.nx+i*grid.nx*grid.ny);
    		}
    for (size_t i=0; i<grid.nx*(grid.ny-1)*grid.nz; i++)
    {
    	kyg.push_back(-kyg[i]);
    }
    //***generate d/dz
    //*Entries (l,j,k)
    for (size_t i=0; i<grid.nz-1; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			lzg.push_back(k+j*grid.nx+i*grid.nx*grid.ny);
    		}
    for (size_t i=0; i<grid.nz-1; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			jzg.push_back(k+j*grid.nx+i*grid.nx*grid.ny);
    		}
    for (size_t i=0; i<grid.nz-1; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			kzg.push_back(-2.0/(1.0*1.0*cellwidth_z[i]+1.0*1.0*cellwidth_z[i+1]));
    		}
    //*Entries (l+1,j,k)
    for (size_t i=0; i<grid.nx*grid.ny*(grid.nz-1); i++)
    {
    	lzg.push_back(lzg[i]);
    }
    for (size_t i=1; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			jzg.push_back(k+j*grid.nx+i*grid.nx*grid.ny);
    		}
    for (size_t i=0; i<grid.nx*grid.ny*(grid.nz-1); i++)
    {
    	kzg.push_back(-kzg[i]);
    }

	//*******Generate index and value of non-zero elements of the gradient matrix S(c) which is a sparse matrix.
	std::vector<int> lxs,lys,lzs,jxs,jys,jzs; //the index of non-zero elements in sparse matrix S.
    std::vector<float> kxs,kys,kzs;           //the value of non-zero elements in sparse matrix S.

	 /*****Diagonal matrix S[nax+nay+naz][nax+nay+naz], nax=(grid.nx-1)*grid.ny*grid.nz,
	 * nay=grid.nx*(grid.ny-1)*grid.nz, naz=grid.nx*grid.ny*(grid.nz-1).
	 *for (size_t i=0; i<nax; i++), S[lxs[i]][jxs[i]]=kxs[i];
	 *for (size_t i=0; i<nay; i++), S[lys[i]+nax][jys[i]+nax]=kys[i]
	 *for (size_t i=0; i<naz; i++), S[lzs[i]+nax+nay][jzs[i]+nax+nay]=kzs[i]
	 *****/

    std::vector<float> rhof_x, rhof_y, rhof_z;
    //***generate x coefficients
    //*Average rho on x face
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=1; k<grid.nx; k++)
    		{
    			rhof_x.push_back((cellwidth_x[k]*cellwidth_y[j]*cellwidth_z[i]*grid.rho[k+j*grid.nx+i*grid.nx*grid.ny]
    			                +cellwidth_x[k-1]*cellwidth_y[j]*cellwidth_z[i]*grid.rho[k-1+j*grid.nx+i*grid.nx*grid.ny])
    					         /(cellwidth_x[k]*cellwidth_y[j]*cellwidth_z[i]+cellwidth_x[k-1]*cellwidth_y[j]*cellwidth_z[i]));
    		}
    for (size_t i=0; i<(grid.nx-1)*grid.ny*grid.nz; i++)
    {
    	lxs.push_back(i);
    	jxs.push_back(i);
    	kxs.push_back(rhof_x[i]);
    }
    //***generate y coefficients
    //*Average rho on y face
    for (size_t i=0; i<grid.nz; i++)
    	for (size_t j=1; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			rhof_y.push_back((cellwidth_x[k]*cellwidth_y[j]*cellwidth_z[i]*grid.rho[k+j*grid.nx+i*grid.nx*grid.ny]
    			                +cellwidth_x[k]*cellwidth_y[j-1]*cellwidth_z[i]*grid.rho[k+(j-1)*grid.nx+i*grid.nx*grid.ny])
    					         /(cellwidth_x[k]*cellwidth_y[j]*cellwidth_z[i]+cellwidth_x[k]*cellwidth_y[j-1]*cellwidth_z[i]));
    		}
    for (size_t i=0; i<grid.nx*(grid.ny-1)*grid.nz; i++)
    {
    	lys.push_back(i);
    	jys.push_back(i);
    	kys.push_back(rhof_y[i]);
    }
    //***generate z coefficients
    //*Average rho on z face
    for (size_t i=1; i<grid.nz; i++)
    	for (size_t j=0; j<grid.ny; j++)
    		for (size_t k=0; k<grid.nx; k++)
    		{
    			rhof_z.push_back((cellwidth_x[k]*cellwidth_y[j]*cellwidth_z[i]*grid.rho[k+j*grid.nx+i*grid.nx*grid.ny]
    			                +cellwidth_x[k]*cellwidth_y[j]*cellwidth_z[i-1]*grid.rho[k+j*grid.nx+(i-1)*grid.nx*grid.ny])
    					         /(cellwidth_x[k]*cellwidth_y[j]*cellwidth_z[i]+cellwidth_x[k]*cellwidth_y[j]*cellwidth_z[i-1]));
    		}
    for (size_t i=0; i<grid.nx*grid.ny*(grid.nz-1); i++)
    {
    	lzs.push_back(i);
    	jzs.push_back(i);
    	kzs.push_back(rhof_z[i]);
    }

    //*******Generate Forward operator A=D*S*G.
    size_t np=grid.nx*grid.ny*grid.nz;
    size_t nax=(grid.nx-1)*grid.ny*grid.nz;
    size_t nay=grid.nx*(grid.ny-1)*grid.nz;
    size_t naz=grid.nx*grid.ny*(grid.nz-1);

    //generate sparse matrix D,
    typedef Eigen::Triplet<float> DT;
    std::vector<DT> DV;
    DV.reserve(2*(nax+nay+naz));
    for (size_t i=0; i<2*nax; i++)
    {
    	DV.push_back(DT(lxd[i],jxd[i],kxd[i]));
    }
    for (size_t i=0; i<2*nay; i++)
    {
    	DV.push_back(DT(lyd[i],jyd[i]+nax,kyd[i]));
    }
    for (size_t i=0; i<2*naz; i++)
    {
    	DV.push_back(DT(lzd[i],jzd[i]+nax+nay,kzd[i]));
    }
    Eigen::SparseMatrix<float> D(np, nax+nay+naz);
    D.setFromTriplets(DV.begin(), DV.end());

    //generate sparse matrix S,
    typedef Eigen::Triplet<float> ST;
    std::vector<ST> SV;
    SV.reserve(nax+nay+naz);
    for (size_t i=0; i<nax; i++)
    {
    	SV.push_back(ST(lxs[i],jxs[i],1.0/kxs[i]));
    }
    for (size_t i=0; i<nay; i++)
    {
    	SV.push_back(ST(lys[i]+nax,jys[i]+nax,1.0/kys[i]));
    }
    for (size_t i=0; i<naz; i++)
    {
    	SV.push_back(ST(lzs[i]+nax+nay,jzs[i]+nax+nay,1.0/kzs[i]));
    }
    Eigen::SparseMatrix<float> S(nax+nay+naz, nax+nay+naz);
    S.setFromTriplets(SV.begin(), SV.end());

    //generate sparse matrix G,
    typedef Eigen::Triplet<float> GT;
    std::vector<GT> GV;
    GV.reserve(2*(nax+nay+naz));
    for (size_t i=0; i<2*nax; i++)
    {
    	GV.push_back(GT(lxg[i],jxg[i],kxg[i]));
    }
    for (size_t i=0; i<2*nay; i++)
    {
    	GV.push_back(GT(lyg[i]+nax,jyg[i],kyg[i]));
    }
    for (size_t i=0; i<2*naz; i++)
    {
    	GV.push_back(GT(lzg[i]+nax+nay,jzg[i],kzg[i]));
    }
    Eigen::SparseMatrix<float> G(nax+nay+naz, np);
    G.setFromTriplets(GV.begin(), GV.end());

    //generate sparse matrix, forward operator A
    Eigen::SparseMatrix<float> A1(np, nax+nay+naz),A(np, np);
    A1=D*S; A=A1*G;
    A.coeffRef(0,0)=1.0/(cellwidth_x[0]*cellwidth_y[0]*cellwidth_z[0]);

    //*****Generate source terms q and receiver terms Q*********************************//
    /* The coordinate of all independent source and receiver pairs for each independent source pair are needed.
     * Solve A*u=q using BICGSTAB (build-in class in Eigen, bi conjugate gradient stabilized solver with preconditioning) to get potential
     * everywhere in the volume for each source pair. Then loop for all source pairs to compute forward data. For each source pair,
     * using Q multiply u (Q*u, Q is terms for all receiver pairs of this source pair) to get forward data.
     * For BICGSTAB, incomplete LU decomposition of A is needed to calculate faster.(if not, is it ok????)
     */
    //generate source term q,
    std::vector<float> q;
    std::vector<float> zpos(grid.nz+1, 0),xpos(grid.nx+1, 0),ypos(grid.ny+1, 0);
    std::vector<float> centerzeroxpos(grid.nx+1, 0),centerzeroypos(grid.ny+1, 0);
    std::vector<float> cellcenterzpos(grid.nz, 0),cellcenterxpos(grid.nx, 0),cellcenterypos(grid.ny, 0);
    //build the 3d grid - numbered from 0 to maximum extent,
    for (size_t i=0; i<grid.nz; i++)
    {
    	zpos[i+1]=zpos[i]+cellwidth_z[i];
    }
    for (size_t i=0; i<grid.nx; i++)
    {
    	xpos[i+1]=xpos[i]+cellwidth_x[i];
    }
    for (size_t i=0; i<grid.ny; i++)
    {
    	ypos[i+1]=ypos[i]+cellwidth_y[i];
    }
    //center the grid about zero. x and y's zero is in the center of 3d volume, z zero at surface.
    for (size_t i=0; i<grid.nx+1; i++)
    {
    	centerzeroxpos[i]=xpos[i]-xpos[grid.nx]/2.0;
    }
    for (size_t i=0; i<grid.ny; i++)
    {
    	centerzeroypos[i]=ypos[i]-ypos[grid.ny]/2.0;
    }
    //set cellcenter's position,
    for (size_t i=0; i<grid.nx; i++)
    {
    	cellcenterxpos[i]=centerzeroxpos[i]+cellwidth_x[i]/2.0;
    }
    for (size_t i=0; i<grid.ny; i++)
    {
    	cellcenterypos[i]=centerzeroypos[i]+cellwidth_y[i]/2.0;
    }
    for (size_t i=0; i<grid.nz; i++)
    {
    	cellcenterzpos[i]=zpos[i]+cellwidth_z[i]/2.0;
    }
    //creat q,
    for (size_t i=0; i<geo.nsource; i++)
    {
    	float sz1_n=geo.sz1[i];
    	float sz2_n=geo.sz2[i];
    	std::vector<float> q1(Interpolate_N(cellcenterxpos,cellcenterypos,cellcenterzpos, geo.sx1[i], geo.sy1[i],sz1_n));
    	std::vector<float> q2(Interpolate_N(cellcenterxpos,cellcenterypos,cellcenterzpos, geo.sx2[i], geo.sy2[i],sz2_n));
    	for (size_t j=0; j<grid.nz; j++)
    		for (size_t k=0; k<grid.ny; k++)
    			for (size_t l=0; l<grid.nz; l++)
    			{
    				q.push_back((q2[l+k*grid.nx+j*grid.nx*grid.ny]-q1[l+k*grid.nx+j*grid.nx*grid.ny])/(cellwidth_x[l]*cellwidth_y[k]*cellwidth_z[j]));
    			}
    	q1.clear();
    	q2.clear();
    }
    //generate receiver term Q, only consider dipole configuration. If not dipole, how to write code??
    std::vector<float> Q;
    for (size_t i=0; i<geo.nsource*geo.nreceiver; i++)
    {
    	std::vector<float> Q1(Interpolate_N(cellcenterxpos,cellcenterypos,cellcenterzpos, geo.rx1[i], geo.ry1[i],geo.rz1[i]));
    	std::vector<float> Q2(Interpolate_N(cellcenterxpos,cellcenterypos,cellcenterzpos, geo.rx2[i], geo.ry2[i],geo.rz2[i]));
    	for (size_t j=0; j<grid.nz; j++)
    		for (size_t k=0; k<grid.ny; k++)
    			for (size_t l=0; l<grid.nz; l++)
    			{
    				Q.push_back(Q2[l+k*grid.nx+j*grid.nx*grid.ny]-Q1[l+k*grid.nx+j*grid.nx*grid.ny]);
    			}
    	Q1.clear();
    	Q2.clear();
    }
    //calculate forward data for every source pair,
    Eigen::VectorXf eu,b,matdata;
    Eigen::MatrixXf eQ;
    std::vector<float> dt(geo.nreceiver, 0.0);
    for (size_t i=0; i<geo.nsource; i++)
    {
        Eigen::BiCGSTAB<SparseMatrix<float>, Eigen::IncompleteLUT<float> >  BCGST;
        BCGST.preconditioner().setDroptol(0.00001);
        BCGST.compute(A);
    	for (size_t j=0; j<np; j++)
    	{
    		b(j)=q[j+i*np];
    	}
    	eu=BCGST.solve(b);
    	for (size_t k=0; k<geo.nreceiver; k++)
    	{
    		for (size_t l=0; l<np; l++)
    		{
    			eQ(k,l)=Q[l+k*np+i*np*geo.nreceiver];
    		}
    	}
    	matdata=eQ*eu;
    	for (size_t m=0; m<geo.nreceiver; m++)
    	{
    		dt[m]=matdata(m);
    		data->dcal.push_back(dt[m]);
    	}
    	dt.clear();
    }

    return(1);
  }
}































