/*
 * time_grav_term.cpp
 *
 *  Created on: Sep 5, 2008
 *      Author: mmoorkamp
 */
#include "../Gravity/ThreeDGravityModel.h"
#include "../Gravity/BasicGravElements.h"
#include <iostream>
#include <boost/date_time/posix_time/posix_time.hpp>
/* Compute gravity effect of a single model element block.
 * Uses the "moment method" that expands the gravitational potential in
 * terms of Legendre polynomials, truncates the series at l=2, sets the
 * element center of mass (COM) at the origin, and then takes the
 * derivative of the potential with respect to z to find the vertical
 * gravitational force.  This has been specialized to a rectangular
 * block with side lengths 2a, 2b, and 2c.
 *
 * Note that for a cube, B02 and B22 are 0 ==> can treat a cube like a
 * sphere and compress the mass to the COM and use Newton's law!
 *
 * Taken from:  Interpretation Theory in Applied Geophysics, by Grant &
 * West. 1965 - pages 222-227
 */

#define GAMMA   6.671e-11

double block(double x, double y, double z, double a, double b, double c,
    double rho)
  {
    double B00, B02, B22;
    double r, r2, r3, r5, r7;
    double F;

    /*
     This computation wants z positive up, but
     rest of model system has z positive down!
     */
    z = -z;

    r = sqrt(x * x + y * y + z * z);
    if (r == 0)
      return (0.0); /* observation point at COM of element
       so ignore this element */
    r2 = r * r;
    r3 = r2 * r;
    r5 = r3 * r2;
    r7 = r5 * r2;

    B00 = 8 * GAMMA * rho * a * b * c;
    B02 = B00 * (2 * c * c - a * a - b * b) / 6.0;
    B22 = B00 * (a * a - b * b) / 24.0;

    F = B00 * z / r3;
    //F += 3* B02 * z * (2 * r2 * r2 - 5 * (x*x - y*y))/(2 * r7);
    F += (B02 / (2 * r5)) * ((5 * z * (3 * z * z - r2) / r2) - 4 * z);
    F += B22 * 15 * z * (x * x - y * y) / r7;
    return F;
  }

int main()
  {
    srand(time(NULL));
    const unsigned int nruns = 1e5;
    double otherbox = 0.0;
    double mybox = 0.0;
    boost::posix_time::ptime mystarttime =
        boost::posix_time::microsec_clock::local_time();
    const double x_meas = rand();
    const double y_meas = rand();
    const double z_meas = -rand();
    const double x_corn = rand();
    const double y_corn = rand();
    const double z_corn = rand();
    const double x_length = rand();
    const double y_length = rand();
    const double z_length = rand();
    for (unsigned int i = 0; i < nruns; ++i)
      {
        mybox = jif3D::CalcGravBoxTerm(x_meas, y_meas, z_meas, x_corn, y_corn,
            z_corn, x_length, y_length, z_length);
      }
    boost::posix_time::ptime myendtime =
        boost::posix_time::microsec_clock::local_time();

    boost::posix_time::ptime newstarttime =
        boost::posix_time::microsec_clock::local_time();
    for (unsigned int i = 0; i < nruns; ++i)
      {
        otherbox = block(x_meas - x_corn - x_length / 2, y_meas - y_corn
            - y_length / 2, z_meas - z_corn - z_length / 2, x_length / 2.0,
            y_length / 2.0, z_length / 2.0, 1) * 1000.0;
      }
    boost::posix_time::ptime newendtime =
        boost::posix_time::microsec_clock::local_time();
    std::cout << mybox << " " << otherbox << " " << fabs(mybox-otherbox)/mybox << std::endl;
    std::cout << "Times: " << myendtime - mystarttime << " " << newendtime
        - newstarttime << std::endl;
    std::cout << "Relative: " << (myendtime - mystarttime).total_microseconds()
        / (newendtime - newstarttime).total_microseconds() << std::endl;
  }

