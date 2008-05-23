	
/*
 * GMT_delaunay performs a Delaunay triangulation on the input data
 * and returns a list of indeces of the points for each triangle
 * found.  Algorithm translated from
 * Watson, D. F., ACORD: Automatic contouring of raw data,
 *   Computers & Geosciences, 8, 97-101, 1982.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "inv3d.h"
#include "fkt.h"

#define FALSE 0
#define TRUE 1

int delaunay (double *x_in, double *y_in, int n, int **link)
              	/* Input point x coordinates */
              	/* Input point y coordinates */
      			/* Number of input points */
            	/* pointer to List of points per triangle.  Vertices for triangle no i is in link[i*3], link[i*3+1], link[i*3+2] */

				/* Output: Number of triangles*/
{
	int ix[3], iy[3];
	int i, j, ij, nuc, jt, km, id, isp, l1, l2, k, k1, jz, i2, kmt, kt, done, size;
	int *index, *stack, *x_tmp, *y_tmp;
	double det[2][3], *x_circum, *y_circum, *r2_circum, *x, *y;
	double xmin, xmax, ymin, ymax, datax, dx, dy, dsq, dd;
	
	size = 10 * n + 1;
	n += 3;
	
	index = (int *) memory (NULL, (3 * size), sizeof (int), "Delaunay");
	stack = (int *) memory (NULL, size, sizeof (int), "Delaunay");
	x_tmp = (int *) memory (NULL, size, sizeof (int), "Delaunay");
	y_tmp = (int *) memory (NULL, size, sizeof (int), "Delaunay");
	x_circum = (double *) memory (NULL, size, sizeof (double), "Delaunay");
	y_circum = (double *) memory (NULL, size, sizeof (double), "Delaunay");
	r2_circum = (double *) memory (NULL, size, sizeof (double), "Delaunay");
	x = (double *) memory (NULL, n, sizeof (double), "Delaunay");
	y = (double *) memory (NULL, n, sizeof (double), "Delaunay");
	
	x[0] = x[1] = -1.0;	x[2] = 5.0;
	y[0] = y[2] = -1.0;	y[1] = 5.0;
	x_circum[0] = y_circum[0] = 2.0;	r2_circum[0] = 18.0;
	
	ix[0] = ix[1] = 0;	ix[2] = 1;
	iy[0] = 1;	iy[1] = iy[2] = 2;
	
	for (i = 0; i < 3; i++) index[i] = i;
	for (i = 0; i < size; i++) stack[i] = i;
	
	xmin = ymin = 1.0e100;	xmax = ymax = -1.0e100;
	
	for (i = 3, j = 0; i < n; i++, j++) {	/* Copy data and get extremas */
		x[i] = x_in[j];
		y[i] = y_in[j];
		if (x[i] > xmax) xmax = x[i];
		if (x[i] < xmin) xmin = x[i];
		if (y[i] > ymax) ymax = y[i];
		if (y[i] < ymin) ymin = y[i];
	}
	
	/* Normalize data */
	
	/*NEW*/
	if((xmax - xmin) >= (ymax - ymin))
		datax = 1.0/(xmax - xmin);
	else
		datax = 1.0/(ymax - ymin);
/*	datax = 1.0 / MAX (xmax - xmin, ymax - ymin);*/
	/*END NEW*/	

	for (i = 3; i < n; i++) {
		x[i] = (x[i] - xmin) * datax;
		y[i] = (y[i] - ymin) * datax;
	}
	
	isp = id = 1;
	for (nuc = 3; nuc < n; nuc++) {
	
		km = 0;
		
		for (jt = 0; jt < isp; jt++) {	/* loop through established 3-tuples */
		
			ij = 3 * jt;
			
			/* Test if new data is within the jt circumcircle */
			
			dx = x[nuc] - x_circum[jt];
			if ((dsq = r2_circum[jt] - dx * dx) < 0.0) continue;
			dy = y[nuc] - y_circum[jt];
			if ((dsq -= dy * dy) < 0.0) continue;

			/* Delete this 3-tuple but save its edges */
			
			id--;
			stack[id] = jt;
			
			/* Add edges to x/y_tmp but delete if already listed */
			
			for (i = 0; i < 3; i++) {
				l1 = ix[i];
				l2 = iy[i];
				if (km > 0) {
					kmt = km;
					for (j = 0, done = FALSE; !done && j < kmt; j++) {
						if (index[ij+l1] != x_tmp[j]) continue;
						if (index[ij+l2] != y_tmp[j]) continue;
						km--;
						if (j >= km) {
							done = TRUE;
							continue;
						}
						for (k = j; k < km; k++) {
							k1 = k + 1;
							x_tmp[k] = x_tmp[k1];
							y_tmp[k] = y_tmp[k1];
							done = TRUE;
						}
					}
				}
				else
					done = FALSE;
				if (!done) {
					x_tmp[km] = index[ij+l1];
					y_tmp[km] = index[ij+l2];
					km++;
				}
			}
		}
			
		/* Form new 3-tuples */
			
		for (i = 0; i < km; i++) {
			kt = stack[id];
			ij = 3 * kt;
			id++;
				
			/* Calculate the circumcircle center and radius
			 * squared of points ktetr[i,*] and place in tetr[kt,*] */
				
			for (jz = 0; jz < 2; jz++) {
				i2 = (jz == 0) ? x_tmp[i] : y_tmp[i];
				det[jz][0] = x[i2] - x[nuc];
				det[jz][1] = y[i2] - y[nuc];
				det[jz][2] = 0.5 * (det[jz][0] * (x[i2] + x[nuc]) + det[jz][1] * (y[i2] + y[nuc]));
			}
			dd = 1.0 / (det[0][0] * det[1][1] - det[0][1] * det[1][0]);
			x_circum[kt] = (det[0][2] * det[1][1] - det[1][2] * det[0][1]) * dd;
			y_circum[kt] = (det[0][0] * det[1][2] - det[1][0] * det[0][2]) * dd;
			dx = x[nuc] - x_circum[kt];
			dy = y[nuc] - y_circum[kt];
			r2_circum[kt] = dx * dx + dy * dy;
			index[ij++] = x_tmp[i];
			index[ij++] = y_tmp[i];
			index[ij] = nuc;
		}
		isp += 2;
	}
	
	for (jt = i = 0; jt < isp; jt++) {
		ij = 3 * jt;
		if (index[ij] < 3 || r2_circum[jt] > 1.0) continue;
		index[i++] = index[ij++] - 3;
		index[i++] = index[ij++] - 3;
		index[i++] = index[ij++] - 3;
	}
	
	index = (int *) memory ((char *)index, i, sizeof (int), "Delaunay");
	*link = index;
	

	free(x_tmp);
	free(y_tmp);
	free(x_circum);
	free(y_circum);
	free(r2_circum);
	free(x);
	free(y);
	free(stack);
	
	return (i/3);
}

#undef FALSE
#undef TRUE

/*-------------------------------------------------------------*/
/* Determine the height Z of the points (X,Y) by means of the known coordinates (x,y,z); The Delaunay triangulation HAD TO  */
/* be performed before and the triangles calculated by the Delaunay are used now for the height interpolation. To assign the right triangle */
/* and to calculate the heights "barycentric coordinates" were used (see in: Contouring (1992) from D.F. Watson; Pergamon Press; pages 76-78)*/
/* Parameter:   *X,*Y := Pointer on the coordinates for which the heights Z should be calculated*/
/*                 *Z := Pointer on the heights Z, which will be calculted by this routine */
/*             numXYZ := Number of the (X,Y) coordinate pairs */
/*           *x,*y,*z := pointer on the coordinates used to build up the triangulation pattern*/
/*             numxyz := Number of the (x,y) coordinate pairs*/
/*   *corner_triangle := pointer on the corner of the triangles; the values i of the vector correspond to the index of x(i), y(i) and z(i) coordinates; */
/*						 all triangles are listed consecutively leading to (3*ntriangles) values */ 
/*        ntriangles  := Number of triangles */

int TriangleInterpol(double *X,double *Y, double *Z, long numXYZ ,double *x ,double *y ,double *z ,int numxyz ,int *corner_triangles,int ntriangles)
{
	long i;
	int j, start_triangle;			/*Triangle, which will be first investigated*/
	double ax,ay,az,bx,by,bz,cx,cy,cz;	/*coordinates of the triangle corners a,b,c */
	double wa,wb,wc;					/*barycentric coordinates*/
	double det;							/*area of the triangle for normalization*/


	start_triangle = 0;

	/*Loop over (X,Y) coordinates*/
	for(i=0;i<numXYZ;i++)
	{

		/*Loop I over triangles*/
		for(j=start_triangle;j<ntriangles;j++)
		{
			/* Calculating the coordinates of the triangle corners*/
			if(corner_triangles[(3*j)+2] > numxyz)
			{
				printf("The positions of the (x,y,z) coordinates do not fit with the triangles\n!");
				exit(0);
			}

			ax = x[corner_triangles[3*j]];
			ay = y[corner_triangles[3*j]];
			az = z[corner_triangles[3*j]];

			bx = x[corner_triangles[(3*j)+1]];
			by = y[corner_triangles[(3*j)+1]];
			bz = z[corner_triangles[(3*j)+1]];

			cx = x[corner_triangles[(3*j)+2]];
			cy = y[corner_triangles[(3*j)+2]];
			cz = z[corner_triangles[(3*j)+2]];

			/* Calculating the barycentric coordinates */
			det=((bx-ax)*(cy-ay))-((cx-ax)*(by-ay)); 
			wa = (((bx - X[i])*(cy - Y[i])) - ((cx - X[i])*(by - Y[i])))/det;
			wb = (((cx - X[i])*(ay - Y[i])) - ((ax - X[i])*(cy - Y[i])))/det;
			wc = (((ax - X[i])*(by - Y[i])) - ((bx - X[i])*(ay - Y[i])))/det;

			/*Is the point in THIS triangle ???*/
			if(wa >= 0 && wb >= 0 && wc >= 0)
			{
				Z[i] = wa*az + wb*bz + wc*cz;
				start_triangle = j;
				goto next;
			}
		}

		/*Loop II over triangles*/
		for(j=0;j<start_triangle;j++)
		{
			/* Calculating the coordinates of the triangle corners*/
			if(corner_triangles[(3*j)+2] > numxyz)
			{
				printf("The positions of the (x,y,z) coordinates do not fit with the triangles\n!");
				exit(0);
			}

			ax = x[corner_triangles[3*j]];
			ay = y[corner_triangles[3*j]];
			az = z[corner_triangles[3*j]];

			bx = x[corner_triangles[(3*j)+1]];
			by = y[corner_triangles[(3*j)+1]];
			bz = z[corner_triangles[(3*j)+1]];

			cx = x[corner_triangles[(3*j)+2]];
			cy = y[corner_triangles[(3*j)+2]];
			cz = z[corner_triangles[(3*j)+2]];
			
			/* Calculating the barycentric coordinates */ 
			det=((bx-ax)*(cy-ay))-((cx-ax)*(by-ay)); 
			wa = (((bx - X[i])*(cy - Y[i])) - ((cx - X[i])*(by - Y[i])))/det;
			wb = (((cx - X[i])*(ay - Y[i])) - ((ax - X[i])*(cy - Y[i])))/det;
			wc = (((ax - X[i])*(by - Y[i])) - ((bx - X[i])*(ay - Y[i])))/det;

			/*Is the point in THIS triangle ???*/
			if(wa >= 0 && wb >= 0 && wc >= 0)
			{
				Z[i] = wa*az + wb*bz + wc*cz;
				start_triangle = j;
				goto next;
			}
		
		}

		/*The point is outside of all triangles:*/
		printf("For the coordinate (%f,%f) no height was assigned\n",X[i],Y[i]);
		Z[i] = -99999.9;

		next:;
	}

	return(1);
}
