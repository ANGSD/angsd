/*
  mod by thorfinn@binf.ku.dk 7 nov 2011.

  The code looks like fortran 2 c translated program of v 2.1
  http://users.eecs.northwestern.edu/~nocedal/lbfgsb.html

  1) Fixed a bug concerning non-consistency in the optimized parameters (accecing noninitialized values)
  2) Is now thread safe. (removed local static variables)
  3) Added a "data" paramater to the findmax_bfgs, such that global variables can be avoided
  4) modified code such that it is c and c++ compilable.
 */

//these three parameters affect precision of optimization.
//m is not recommended to be higher than 20
//decreasing factr, pgtol will increase precision
//factr is the multiple of machine precision that result will be
//pgtol is size of gradient on exit
#define MVAL 10
#define FACTR 1.0e6
#define PGTOL 1.0e-3


//nbd is a vector of integers of dimension numpars.
//nbd[i]=0 if there are no bounds for parameter i, 
//      =1 if there are only lower bounds
//      =2 if there are both lower/upper bounds
//      =3 if there is only upper bound                
//or send nbd=NULL is equivalent to nbd=2 for all parameters
//
//noisy=0 => no output, noisy=1 => one line of output, noisy < 99 some output,
//noisy>=100 probably too much
//
//dfun is derivative function or send NULL to use numerical derivative
//(getgradient function)

/*
  numpars = dimension of optimization
  invec   = Must have length numpars, will be used at startpoint for optimization. 
            Contains the optimized parameters
  dats    = This can contain a pointer to a generic structure which contains need data
            Can be NULL, if no "global" variables are needed
  fun     = function to optimize
  dfun    = derivative function
  lowbound= lower bounds
  upbound = upper bounds
  nbd,noisy   = see above
  
 */
double findmax_bfgs(int numpars, double *invec,const void*dats, double (*fun)(const double x[],const void*),
		    void (*dfun)(const double x[], double y[]),
		    double *lowbound, double *upbound,
		    int *nbd, int noisy);
/*
void getgradient(int npar, const double invec[],double outvec[],
		 const void*dats,double(*func)(const double [],const void*), 
		 const double* lowbound, const double* upbound);
*/
