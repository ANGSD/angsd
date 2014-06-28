#include <cmath>

double a1f(int nsam)
{
double a1;
int i;
a1 = 0.0;
for (i=1; i<=nsam-1; i++) a1 += 1.0/i;
return (a1);
}


double a2f(int nsam) 
{
double a2;
int i;
a2 = 0.0;
for (i=1; i<=nsam-1; i++) a2 += 1.0/(i*i);
return (a2);
}


double b1f(int nsam){
  double b1;
  b1 = (nsam + 1.0)/(3.0*(nsam-1.0));
  return (b1);
}


double b2f(int nsam){
  double b2;
  b2 = (2*(nsam*nsam + nsam + 3.0))/(9*nsam*(nsam - 1));
  return (b2);
}


double e1f(double a1, double c1) {
  double e1;
  e1 = c1/a1;
  return (e1);
}

double e2f(double a1, double a2, double c2){ 
  double e2;
  e2 = c2/((a1*a1)+a2);
  return (e2);
}


double c1f(double a1, double b1) {
  double c1;
  c1 = b1 - (1/a1);
  return (c1);
}


double c2f(int nsam, double a1, double a2, double b2) {
  double c2;
  c2 = b2 - ((nsam+2)/(a1*nsam)) + (a2/(a1 * a1));
  return (c2);
}




double tajd(int nsam, double thetaW, double sumk){
  static double a1 = a1f(nsam);

  double segsites  = thetaW*a1;
  if( segsites == 0 ) return( 0.0) ;
  
  
  static double a2 = a2f(nsam);
  static double b1 = b1f(nsam);
  static double b2 = b2f(nsam);
  static double c1 = c1f(a1, b1);
  static double c2 = c2f(nsam, a1, a2, b2);
  static double e1 = e1f(a1, c1);
  static double e2 = e2f(a1, a2, c2);

  return( (sumk - (thetaW))/sqrt((e1*segsites) + ((e2*segsites)*(segsites-1))) ) ;
}




double cn(int n){
  double top =  2*n*a1f(n)-4*(n-1);
  double bot =  (n-1)*(n-2);
  return(top/bot);
}



double vd(int n){
  double led1 = a1f(n)*a1f(n)/(a2f(n)+a1f(n)*a1f(n));
  double led2 = cn(n)-(n+1)/(1.0*(n-1));
  return(1+led1*led2);

}


double ud(int n){
  return(a1f(n)-1-vd(n));
  
}


double vf(int n){
  double top = cn(n)+2*(n*n+n+3)/(1.0*(9*n*(n-1)))-2/(1.0*(n-1));
  double bot = a1f(n)*a1f(n)+a2f(n);
  
  return(top/bot);
}

double uf(int n){
  double top = 1+(n+1)/(3.0*(n-1))-4*((n+1)/(1.0*(n-1)*(n-1)))*(a1f(n+1)-(2*n)/(1.0*(n+1)));
  double bot = a1f(n);
  return(top/(bot)-vf(n)  );
}



double fulid(int n, double thetaW, double thetaFL){


  double S = thetaW*a1f(n);
  double L = thetaFL*a1f(n);

  double top = S-L;
  // double bot = ud(n)*S+vd(n)*S*S;
  static double ud_val = ud(n);
  static double vd_val = vd(n);
  double bot = ud_val*S+vd_val*S*S;

  //  fprintf(stderr,"S=%f L=%f top=%f bot=%f\n",S,L,top,bot);
  //  fflush(stderr);
  return(top/sqrt(bot));
}

double fulif(int n, double thetaW, double thetaFL,double thetaPi){

  double S = thetaW*a1f(n);
  double top = thetaPi-thetaFL;
  //  double bot = uf(n)*S+vf(n)*S*S;
  
  static double uf_val = uf(n);
  static double vf_val = vf(n);
  double bot = uf_val*S+vf_val*S*S;
  
  return (top/sqrt(bot));
    
}





double H_var(int n){

  double top = 18*n*n*(3*n+2)*a2f(n+1)-(88*n*n*n+9*n*n-13*n+6);
  double bot = 9*n*(n-1)*(n-1);
  return (top/bot);
}

double fayh (int n, double thetaW, double thetaH,double thetaPi){
  double en = thetaH;
  double to = thetaPi;
  static double a1f_val = a1f(n); 
  double S  =  thetaW*a1f_val;

  double led1 = (n-2)/(6.0*(n-1))*S;
  double led2 = H_var(n)*S*S;
  return ((to-en)/sqrt(led1+led2));

}



double E_var1(int n){
  return(n/(2.0*(n-1))-1/a1f(n));
}

double E_var2(int n){
  double led1 = a2f(n)/(a1f(n)*a1f(n));
  double tmp = n/(1.0*(n-1));
  double led2 = 2*tmp*tmp*a2f(n);
  double led3 = (2*(n*a2f(n)-n+1)/(1.0*(n-1)*a1f(n)));
  double led4 = (3*n+1)/(1.0*(n-1));
  return (led1+led2-led3-led4);
}
  


double zenge(int n, double thetaW, double thetaL){

  double S = thetaW*a1f(n);

  double top =thetaL-thetaW;
  double bot = E_var1(n)*S+E_var2(n)*S*S;
  return(top/sqrt(bot));
}

