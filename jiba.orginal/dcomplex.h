typedef struct DCOMPLEX {long double r,i;} dcomplex;

/*Self-written routines*/
extern dcomplex CBtanh(dcomplex z);
extern int CDerivQuo(double Ar, double Ai, double Br, double Bi, double dAr, double dAi, double dBr, double dBi, double *dCr, double *dCi);

extern dcomplex Cadd(dcomplex a,dcomplex b);
extern dcomplex Csub(dcomplex a,dcomplex b);
extern dcomplex Cmul(dcomplex a,dcomplex b);
extern dcomplex Complex(double re, double im);
extern dcomplex Conjg(dcomplex z);
extern dcomplex Cdiv(dcomplex a,dcomplex b);
extern long double Cabs(dcomplex z);
extern long double Cabs_sq(dcomplex z);
extern dcomplex Csqrt(dcomplex z);
extern dcomplex RCmul(double x,dcomplex a);
extern dcomplex rect2polar(dcomplex z);
extern dcomplex polar2rect(dcomplex polar);
extern dcomplex Cpow(dcomplex z,int n);
extern dcomplex Cexp(dcomplex z);
extern dcomplex Ccexp(double a, double b);

