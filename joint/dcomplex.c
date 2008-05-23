/* d_complex.c 

   Author: Numerical recipes
   Revised by: Jon Hamkins
*/

#include <math.h>
#include "dcomplex.h"

dcomplex Cadd(dcomplex a,dcomplex b)
{   dcomplex c;
    c.r=a.r+b.r;
    c.i=a.i+b.i;
    return c;
}

dcomplex Csub(dcomplex a,dcomplex b)
{   dcomplex c;
    c.r=a.r-b.r;
    c.i=a.i-b.i;
    return c;
}

dcomplex Cmul(dcomplex a,dcomplex b)
{   dcomplex c;
    c.r=a.r*b.r-a.i*b.i;
    c.i=a.i*b.r+a.r*b.i;
    return c;
}

dcomplex Complex(double re,double im)
{   dcomplex c;
    c.r=re;
    c.i=im;
    return c;
}

dcomplex Conjg(dcomplex z)
{   dcomplex c;
    c.r=z.r;
    c.i = -z.i;
    return c;
}

dcomplex Cdiv(dcomplex a,dcomplex b)
{   dcomplex c;
    long double r,den;
    if (fabsl(b.r) >= fabsl(b.i)) {
        r=b.i/b.r;
        den=b.r+r*b.i;
        c.r=(a.r+r*a.i)/den;
        c.i=(a.i-r*a.r)/den;
    } else {
        r=b.r/b.i;
        den=b.i+r*b.r;
        c.r=(a.r*r+a.i)/den;
        c.i=(a.i*r-a.r)/den;
    }
    return c;
}

long double Cabs(dcomplex z)     /* magnitude of complex number */
                            /* compute in manner that avoids overflow errors */
{   long double x,y,ans,temp;
    x=fabsl(z.r);
    y=fabsl(z.i);
    if (x == 0.0)
        ans=y;
    else if (y == 0.0)
        ans=x;
    else if (x > y) {
        temp=y/x;
        ans=x*sqrt(1.0+temp*temp);
    } else {
        temp=x/y;
        ans=y*sqrt(1.0+temp*temp);
    }
    return ans;
}

long double Cabs_sq(dcomplex z)
{   long double x,y;
    x=z.r;
    y=z.i;
    return x*x+y*y;
}

dcomplex Csqrt(dcomplex z)
{   dcomplex c;
    long double x,y,w,r;
    if ((z.r == 0.0) && (z.i == 0.0)) {
        c.r=0.0;
        c.i=0.0;
        return c;
    } else {
        x=fabsl(z.r);
        y=fabsl(z.i);
        if (x >= y) {
            r=y/x;
            w=sqrtl(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
        } else {
            r=x/y;
            w=sqrtl(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
        }
        if (z.r >= 0.0) {
            c.r=w;
            c.i=z.i/(2.0*w);
        } else {
            c.i=(z.i >= 0) ? w : -w;
            c.r=z.i/(2.0*c.i);
        }
        return c;
    }
}

dcomplex RCmul(double x,dcomplex a)
{   dcomplex c;
    c.r=x*a.r;
    c.i=x*a.i;
    return c;
}

dcomplex rect2polar(dcomplex z)
/*************************************************************************
*   Converts complex numbers from rectangular coordinates to polar
*
*   INPUTS:  z in rect format
*   OUTPUTS: none
*   CALLS:   none
*   RETURNS: polar form of z, with |z| in the .r and <z in the .i
*
**************************************************************************/
{
dcomplex polar;

polar.i=atan2l(z.i,z.r);
polar.r=sqrtl(z.r*z.r+z.i*z.i);

return polar;
}

dcomplex polar2rect(dcomplex polar)
/*************************************************************************
*   Converts complex numbers from polar coordinates to rectangular
*
*   INPUTS:  z in polar format, with |z| in the .r and <z in the .i
*   OUTPUTS: none
*   CALLS:   none
*   RETURNS: z in standard rectangular format.
*
**************************************************************************/
{
int negative=0;
long double temp;
dcomplex rect;

#define TWO_PI 6.28318530718
#define M_PI 3.14159265358979

temp=fmodl(polar.i,TWO_PI);
if(temp<0.0)
    {
    if(temp>-M_PI)
        negative=1;
    }
else
    {
    if(temp>M_PI)
        negative=1;
    }
temp=cosl(temp);
rect.r=polar.r*temp;
rect.i=polar.r*sqrtl(1.0-temp*temp);
if(negative) rect.i*=-1.0;

return rect;
}

dcomplex Cpow(dcomplex z,int n)
/*************************************************************************
*   Computes integer powers of complex numbers
*
*   ALGORITHM:  convert to polars, find z^n, convert back to rectangular
*   INPUTS:  z and n
*   OUTPUTS: none
*   CALLS:   polar2rect, rect2polar()
*   RETURNS: z^n
*   NOTE: for small vales of n, it may be better to compute using the
*         binomial expansion - avoids trig calls
*
**************************************************************************/
{
dcomplex    z2n;
if(n==1) return z;
if(!n) {z2n.r=1.0; z2n.i=0.0; return z2n;}

z2n=rect2polar(z);
if(z2n.r==0.0) {z2n.r=z2n.i=0.0; return z2n;}
z2n.r=powl(z2n.r,(double) n);
z2n.i=z2n.i*n;
z2n=polar2rect(z2n);
return z2n;
}

/*************************************************************************
*
* Cexp
*
* Computes e^z, where z is a complex number
* 
* Method: e^z = e^(a+ib) = e^a e^(ib) 
              = e^a ( cos(b) + i sin(b))
              = e^a cos(b) + i e^a sin(b)
**************************************************************************/
dcomplex Cexp(dcomplex z)
{
	dcomplex ans;
	register long double e;

	e = expl(z.r);           /* compute only once to save time */
	ans.r = e * cosl(z.i);   /* real part of e^z */
	ans.i = e * sinl(z.i);   /* imaginary part of e^z */
	return(ans);
}

/*************************************************************************
*
* Ccexp
*
* Computes a e^(ib), where a and b are real numbers
* 
* Method: a e^(ib) = a ( cos(b) + i sin(b)) = a cos(b) + i a sin(b)
**************************************************************************/
dcomplex Ccexp(double a, double b)
{
	dcomplex ans;

	ans.r = a * cos(b);   /* real part of e^z */
	ans.i = a * sin(b);   /* imaginary part of e^z */
	return(ans);
}
/*************************************************************************
*
* CBtanh
*
* Computes the Tangens Hyperbolicus from a complex number z
* (This routine is written from Bjoern Heincke)
**************************************************************************/
dcomplex CBtanh(dcomplex z)
{
	#define LIMIT 600

	dcomplex tmp1,tmp2,zaehler,nenner,ans;
	long double r,den;
	long double e_plus,e_minus;

	/*Make sure that the exponential function can ve calculated*/
	if(z.r > LIMIT)
		z.r = LIMIT;
	if(z.r < ((-1) * LIMIT))
		z.r = (-1) * LIMIT;

	/*Calculate e^z*/
	e_plus = expl(z.r);				/* compute only once to save time */
	tmp1.r = e_plus * cosl(z.i);		/* real part of e^z */
	tmp1.i = e_plus * sinl(z.i);		/* imaginary part of e^z */

	/*Calculate e^(-z)*/
	e_minus = expl(-z.r);			/* compute only once to save time */
	tmp2.r = e_minus * cosl(-z.i);   /* real part of e^-z */
	tmp2.i = e_minus * sinl(-z.i);   /* imaginary part of e^-z */

	/*Calculate the nominator e^z - e^(-z)*/
	zaehler.r = tmp1.r - tmp2.r;
    zaehler.i = tmp1.i - tmp2.i;

	/*Calculate the denominator e^z + e^(-z)*/
	nenner.r = tmp1.r + tmp2.r;
    nenner.i = tmp1.i + tmp2.i;

	/*Divide the nominator by the denominator*/
    if (fabsl(nenner.r) >= fabsl(nenner.i)) 
	{
        r = nenner.i/nenner.r;
        den = nenner.r + r * nenner.i;
        ans.r = (zaehler.r + r*zaehler.i)/den;
        ans.i = (zaehler.i - r*zaehler.r)/den;
    } 
	else 
	{
        r = nenner.r/nenner.i;
        den = nenner.i + r * nenner.r;
        ans.r = (zaehler.r*r+zaehler.i)/den;
        ans.i = (zaehler.i*r-zaehler.r)/den;
    }

	#undef LIMIT

	return(ans);
}

/*************************************************************************
*
* ComplexDerivQuo
*
* Computes the complex derivative dC/dx  with C=A/B and A,B,C are complex and x is real
* -Ar and Br are the real parts of A and B
* -Ai and Bi are the imaginary parts of A and B
* -dAr, dBr and dCr are the real part of (dA/dx), (dB/dx) and (dC/dx) 
* -dAi, dBi and dCi are the imaginary part of (dA/dx), (dB/dx) and (dC/dx) 
* (This routine is written from Bjoern Heincke)
**************************************************************************/
#define EPS 1E-15

extern int CDerivQuo(double Ar, double Ai, double Br, double Bi, double dAr, double dAi, double dBr, double dBi, double *dCr, double *dCi)
{
	double m1,m2,m3,m4,m5,m6;
	double nominator, denominator; 

	/*Calculate the real part of the dC/dx*/
	m1 = (Br*Br + Bi*Bi);
	m2 = ((dAr*Br)+(Ar*dBr)+(dAi*Bi)+(Ai*dBi));
	m3 = (2*Br*dBr+2*Bi*dBi);
	m4 = (Ar*Br + Ai*Bi);

	nominator = m1*m2 - m3*m4;
	denominator = m1*m1;

	if(denominator == 0.0)
		denominator = EPS;

	(*dCr) = (nominator/denominator);

	/*Calculate the imaginary part of the dC/dx*/
	m5 = (double)((dAi*Br)+(Ai*dBr)-(dBi*Ar)-(Bi*dAr));
	m6 = (double)(Ai*Br - Bi*Ar);

	nominator = m1*m5 - m3*m6;

	(*dCi) = (nominator/denominator);

	return(0);
}

#undef EPS