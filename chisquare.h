
#pragma once
#ifndef _chisq_
#define _chisq_

#include <cstdlib>
/*
  taken from numerical recipis 2007 third edition
  modified by thorfinn@binf.ku.dk 30 aug 2013

  ppl on the internet complain that numerecial recp can't be used in academic gnu software, will need to check that and use a different impl.
 */

namespace chisq{
struct Gauleg18 {
  static const int ngau = 18;
  static const double y[18];
  static const double w[18];
};

struct Gamma : Gauleg18 {
  static const int ASWITCH=100;
  static const double EPS;
  static const double FPMIN;
  double gln;
  double gammp(const double a, const double x);
  double gammq(const double a, const double x);
  double gser(const double a, const double x);
  double gcf(const double a, const double x);
  double gammpapprox(double a, double x, int psig);
  double invgammp(double p, double a);
  
};

  double gammln(const double xx);
}

struct Chisqdist : chisq::Gamma {
  double nu,fac;
  Chisqdist(double nnu) : nu(nnu) {
    fac = 0.693147180559945309*(0.5*nu)+chisq::gammln(0.5*nu);
  }
  double p(double x2) {
    if (x2 <= 0.) throw("bad x2 in Chisqdist");
    return exp(-0.5*(x2-(nu-2.)*log(x2))-fac);
  }
  double cdf(double x2) {
    //fprintf(stderr,"x2:%f\n",x2);
    if (x2 < 0.) {
      fprintf(stderr,"bad x2 in Chisqdist:%f \n",x2);
      exit(0);
    }
    return gammp(0.5*nu,0.5*x2);
  }
  double invcdf(double p) {
    
      if (p < 0. || p >= 1.) throw("bad p in Chisqdist");
    return 2.*invgammp(p,0.5*nu);
  }
};

#endif
