
#include <cmath>
#include <iostream>
#include <limits>
#include <cstdio>
#include <cassert>
#include "chisquare.h"

using namespace chisq;
/*
  taken from numerical recipis 2007 third edition
  modified by thorfinn@binf.ku.dk 30 aug 2013

  ppl on the internet complain that numerecial recp can't be used in academic gnu software, will need to check that and use a different impl.
 */


const double Gauleg18::y[18] = {0.0021695375159141994,
				0.011413521097787704,0.027972308950302116,0.051727015600492421,
				0.082502225484340941, 0.12007019910960293,0.16415283300752470,
				0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
				0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
				0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
				0.87126389619061517, 0.95698180152629142};
const double Gauleg18::w[18] = {0.0055657196642445571,
				0.012915947284065419,0.020181515297735382,0.027298621498568734,
				0.034213810770299537,0.040875750923643261,0.047235083490265582,
				0.053244713977759692,0.058860144245324798,0.064039797355015485,
				0.068745323835736408,0.072941885005653087,0.076598410645870640,
				0.079687828912071670,0.082187266704339706,0.084078218979661945,
				0.085346685739338721,0.085983275670394821};

double chisq::gammln(const double xx) {
  int j;
  double x,tmp,y,ser;
  static const double cof[14]={57.1562356658629235,-59.5979603554754912,
			     14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
			       .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
			       -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
			       .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
  if (xx <= 0) throw("bad arg in gammln");
  y=x=xx;
  tmp = x+5.24218750000000000; 
  tmp = (x+0.5)*log(tmp)-tmp;
  ser = 0.999999999999997092;
  for (j=0;j<14;j++) 
    ser += cof[j]/++y;
  return tmp+log(2.5066282746310005*ser/x);
}



double Gamma::gammp(const double a, const double x) {
    if (x < 0.0 || a <= 0.0) 
      throw("bad args in gammp");
    if (x == 0.0) return 0.0;
    else if ((int)a >= ASWITCH) return gammpapprox(a,x,1);
    else if (x < a+1.0) return gser(a,x);
    else return 1.0-gcf(a,x);
}
double Gamma::gammq(const double a, const double x) {

  if (x < 0.0 || a <= 0.0) throw("bad args in gammq");
    if (x == 0.0) return 1.0;
    else if ((int)a >= ASWITCH) 
      return gammpapprox(a,x,0);
    else if (x < a+1.0) return 1.0-gser(a,x);
    else return gcf(a,x); 
}
double Gamma::gser(const double a, const double x) {
    double sum,del,ap;
    gln=gammln(a);
    ap=a;
    del=sum=1.0/a;
    for (;;) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	return sum*exp(-x+a*log(x)-gln);
      }
    }
  }
double Gamma::gcf(const double a, const double x) {
    int i;
    double an,b,c,d,del,h;
    gln=gammln(a);
    b=x+1.0-a;
    c=1.0/FPMIN;
    d=1.0/b;
    h=d;
    for (i=1;;i++) {
      an = -i*(i-a);
      b += 2.0;
      d=an*d+b;
      if (fabs(d) < FPMIN)
	d=FPMIN;
      c=b+an/c;
      if (fabs(c) < FPMIN)
	c=FPMIN;
      d=1.0/d;
      del=d*c;
      h *= del;
      if (fabs(del-1.0) <= EPS)
	break;
    }
    return exp(-x+a*log(x)-gln)*h;
  }
double Gamma::gammpapprox(double a, double x, int psig) {
    int j;
    double xu,t,sum,ans;
    double a1 = a-1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
    gln = gammln(a);
    
    if (x > a1) xu = std::max(a1 + 11.5*sqrta1, x + 6.0*sqrta1);
    else xu = std::max(0.,std::min(a1 - 7.5*sqrta1, x - 5.0*sqrta1));
    sum = 0;
    for (j=0;j<ngau;j++) { 
      t = x + (xu-x)*y[j];
      sum += w[j]*exp(-(t-a1)+a1*(log(t)-lna1));
    }
    ans = sum*(xu-x)*exp(a1*(lna1-1.)-gln);
    return (psig?(ans>0.0? 1.0-ans:-ans):(ans>=0.0? ans:1.0+ans));
  }

double Gamma::invgammp(double p, double a) {
  int j;
  double x,err,t,u,pp,lna1,afac,a1=a-1;
  const double EPS=1.e-8;
  gln=gammln(a);
  if (a <= 0.) throw("a must be pos in invgammap");
  if (p >= 1.) return std::max(100.,a + 100.*sqrt(a));
  if (p <= 0.) return 0.0;
  if (a > 1.) {
      lna1=log(a1);
    afac = exp(a1*(lna1-1.)-gln);
    pp = (p < 0.5)? p : 1. - p;
    t = sqrt(-2.*log(pp));
    x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
    if (p < 0.5) x = -x;
    x = std::max(1.e-3,a*pow(1.-1./(9.*a)-x/(3.*sqrt(a)),3));
  } else {
      t = 1.0 - a*(0.253+a*0.12);
      if (p < t) x = pow(p/t,1./a);
      else x = 1.-log(1.-(p-t)/(1.-t));
  }
  for (j=0;j<12;j++) {
    if (x <= 0.0) return 0.0;
    err = gammp(a,x) - p;
    if (a > 1.) t = afac*exp(-(x-a1)+a1*(log(x)-lna1));
    else t = exp(-x+a1*log(x)-gln);
    u = err/t;
    x -= (t = u/(1.-0.5*std::min(1.,u*((a-1.)/x - 1))));
    if (x <= 0.) x = 0.5*(x + t); 
    if (fabs(t) < EPS*x ) break;
  }
return x;
}



const double Gamma::EPS = std::numeric_limits<double>::epsilon();
const double Gamma::FPMIN = std::numeric_limits<double>::min()/EPS;

/*
//in .h now
struct Chisqdist : Gamma {
  double nu,fac;
  Chisqdist(double nnu) : nu(nnu) {
    fac = 0.693147180559945309*(0.5*nu)+gammln(0.5*nu);
  }
  double p(double x2) {
    if (x2 <= 0.) throw("bad x2 in Chisqdist");
    return exp(-0.5*(x2-(nu-2.)*log(x2))-fac);
  }
  double cdf(double x2) {
    if (x2 < 0.) throw("bad x2 in Chisqdist");
    return gammp(0.5*nu,0.5*x2);
  }
  double invcdf(double p) {
    
      if (p < 0. || p >= 1.) throw("bad p in Chisqdist");
    return 2.*invgammp(p,0.5*nu);
  }
};
*/
#ifdef __WITH_MAIN__
int main(int argc, char **argv){
  fprintf(stderr,"Thorfinn chisquare function\n------------------\n");
  double df;
  std::cout << "Type degrees of freedom"<<std::endl;
  std::cin >> df;
  assert(df>0);
  int type =0;
  
  std::cout << "Type 1-3.\n\t1: density (R::dchisq)\n\t2: cdf (R::pchisq)\n\t3: invcdf (R::qchisq)"<<std::endl;
  std::cin>>type;
  if(type==0) return 0;
  
  Chisqdist tmp(df);

  while(1){
    double val;
    std::cout << "Type point for evaluation"<<std::endl;
    std::cin >>val;
    if(type==1)
      fprintf(stderr,"Chisquare.p(%f,%f)=%e\n",val,df,tmp.p(val));
    if(type==2)
      fprintf(stderr,"Chisquare.cdf(%f,%f)=%e\n",val,df,tmp.cdf(val));
    if(type==3)
      fprintf(stderr,"Chisquare.invcdf(%f,%f)=%e\n",val,df,tmp.invcdf(val));
      
  }

  return 0;
}
#endif
