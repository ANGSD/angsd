void ComputeP11(unsigned numWind,int tk_l,double *P1,double *PP1,double **fw,double **bw,double *stationary){
  for (unsigned i = 0; i < tk_l; i++){
    PP1[i] = 0;
    for (unsigned l = 1; l < numWind; l++)
      PP1[i] += fw[i][l]*P1[i]*bw[l+1][i]/stationary[i];//NOTE: In appendix of K.H. paper it seems to be an extra emission probability for site l+1, it is already inside bw[]
  }
}

void ComputeP22(unsigned numWind,int tk_l,double **P,double **PP,double **fw,double **bw,double *stationary){
  double R1[tk_l];
  double R2[tk_l];

  for (unsigned i = 0; i < tk_l; i++)
    PP[2][i] = 0;
  for (unsigned l = 1; l < numWind; l++){
    

    R1[tk_l - 1] = 0;
    for (int i = tk_l - 2; i >= 0 ; i--)
      R1[i] = R1[i+1] + fw[i+1][l];
    double tmp = 0;
    for (unsigned i = 0; i < tk_l ; i++){
      R2[i] = tmp*P[2][i]+fw[l][i]*P[6][i]+R1[i]*P[7][i];
      tmp = R2[i];
    }
    for (unsigned i = 1; i < tk_l; i++)
      PP[2][i] += R2[i-1]*P[2][i]*bw[i][l+1]/stationary[i];
  }
}

void ComputeP33(unsigned numWind,int tk_l,double *P3,double *PP3,double **fw,double **bw,double *stationary){
  double R1[tk_l];
  for (unsigned i = 0; i < tk_l; i++)
    PP3[i] = 0;
  for (unsigned l = 1; l < numWind; l++){
    R1[tk_l - 1] = 0;
    for (int i = tk_l - 2; i >= 0 ; i--)
      R1[i] = R1[i+1] + fw[i+1][l];
    for (unsigned i = 0; i < tk_l - 1; i++)
      PP3[i] += R1[i]*P3[i]*bw[i][l]/stationary[i];
  }
}

void ComputeP44(unsigned numWind,int tk_l,double *P4,double *PP4,double **fw,double **bw,double *stationary){
  for (unsigned i = 0; i < tk_l; i++){
    PP4[i] = 0;
    for (unsigned l = 1; l < numWind; l++)
      PP4[i] += fw[i][l]*P4[i]*bw[i][l+1]/stationary[i];
  }
}

void ComputeP55(unsigned numWind,int tk_l,double **P,double **PP,double **fw,double **bw,double *stationary){
  double R1[tk_l];
  double R2[tk_l];
  double bR1[tk_l];

  for (unsigned l = 1; l < numWind; l++){
    R1[tk_l - 1] = 0;
    for (int i = tk_l - 2; i >= 0 ; i--)
      R1[i] = R1[i+1] + fw[i+1][l];
  }
  for (unsigned i = 0; i < tk_l; i++)
    PP[5][i] = 0;
  for (unsigned l = 1; l < numWind; l++){
    double tmp = 0;
    for (unsigned i = 0; i < tk_l ; i++){
      R2[i] = tmp*P[2][i]+fw[l][i]*P[6][i]+R1[i]*P[7][i];
      tmp = R2[i];
    }
    bR1[tk_l - 1] = 0;
    for (int i = tk_l - 2; i >= 0 ; i--)
      bR1[i] = bR1[i+1] + bw[i+1][l+1];
    for (unsigned i = 0; i < tk_l - 1; i++)
      PP[5][i] += R2[i]*P[5][i]*bR1[i]/P[0][i];//<- CHECK ptgi
  }
}

void ComputeP66(unsigned numWind,int tk_l,double **P,double **PP,double **fw,double **bw,double *stationary){
  double bR1[tk_l];
  for (unsigned i = 0; i < tk_l; i++)
    PP[6][i] = 0;
  for (unsigned l = 1; l < numWind; l++){
    bR1[tk_l - 1] = 0;
    for (int i = tk_l - 2; i >= 0 ; i--)
      bR1[i] = bR1[i+1] + bw[i+1][l+1];
    for (unsigned i = 0; i < tk_l - 1; i++)
      PP[6][i] += fw[i][l]*P[6][i]*bR1[i]/P[0][i];//<- CHECK btgi
  }
}

void ComputeP77(unsigned numWind,int tk_l,double **P,double **PP,double **fw,double **bw,double *stationary){
  double R1[tk_l];
  double bR1[tk_l];
  for (unsigned i = 0; i < tk_l; i++)
    PP[7][i] = 0;
  for (unsigned l = 1; l < numWind; l++){
    R1[tk_l - 1] = 0;
    for (int i = tk_l - 2; i >= 0 ; i--)
      R1[i] = R1[i+1] + fw[i+1][l];
    
    bR1[tk_l - 1] = 0;
    for (int i = tk_l - 2; i >= 0 ; i--)
      bR1[i] = bR1[i+1] + bw[i+1][l+1];

    for (unsigned i = 0; i < tk_l - 1; i++)
      PP[7][i] += R1[i]*P[7][i]*bR1[i]/P[0][i];//<-CHECK ptgi
  }
}

//TODO: tk_l is in fact the number of time intervals
void ComputeP1(double *tk,int tk_l,double *P,double *epsize,double rho){ 
  for (unsigned i = 0; i < tk_l-1; i++){
    P[i] = 1.0/(1.0+epsize[i]*2.0*rho);
    P[i] *= exp( -rho*2.0*tk[i] ) - exp(-rho*2.0*tk[i+1]-(tk[i+1]-tk[i])/epsize[i]);
    P[i] /= 1.0 - exp( -(tk[i+1]-tk[i])/epsize[i] );
  }
  //Last interval ends with +infinity
  unsigned i = tk_l - 1;
  P[i] = 1.0/(1.0+epsize[i]*2.0*rho)* exp( -rho*2.0*tk[i] );
}

void ComputeP5(double *tk,int tk_l,double *P,double *epsize){
  for (unsigned i = 0; i < tk_l-1; i++)
    P[i] = exp( -(tk[i+1] - tk[i])/epsize[i] );
  P[tk_l-1] = 0.0;
}


void ComputeP6(double *tk,int tk_l,double *P,double *epsize,double rho){
    for (unsigned i = 0; i < tk_l-1; i++){
      P[i] = 1/(1-exp(-(tk[i+1]-tk[i])/epsize[i]));
      P[i] *= exp(-(tk[i+1]-tk[i])/epsize[i]);
      double tmp = exp(-2*rho*tk[i]);
      tmp -= 1/(1-2*rho*epsize[i])*exp(-2*rho*tk[i+1]);
      tmp += 2*rho*epsize[i]/(1 - 2*rho*epsize[i])*exp(-2*rho*tk[i]-(tk[i+1]-tk[i])/epsize[i]);
      P[i] *= tmp;
    }
    P[tk_l - 1] = 0.0;
  }


void ComputeP2(int tk_l,double *P2,double *P5){
    for (unsigned i = 0; i < tk_l; i++)
      P2[i] = 1.0 - P5[i];
  }



void ComputeP3(double *tk,int tk_l,double *P3,double *epsize,double rho){
  for (unsigned i = 0; i < tk_l - 1; i++){
    P3[i] = exp(-tk[i]*2.0*rho);
    P3[i] += epsize[i]*2.0*rho/(1.0 - epsize[i]*2.0*rho)*exp(-(tk[i+1]-tk[i])/epsize[i]-tk[i]*2.0*rho);
    P3[i] -= 1.0/(1.0 - epsize[i]*2.0*rho)*exp(-tk[i+1]*2.0*rho);
  }
  unsigned i = tk_l - 1;
  P3[i] = exp(-tk[i]*2.0*rho);
}


void ComputeP4(double *tk,int tk_l,double *P4,double *epsize,double rho){
  for (unsigned i = 0; i < tk_l-1; i++){
    P4[i] = 1.0/(1.0 - exp(-(tk[i+1]-tk[i])/epsize[i]) );
    double tmp = 2.0*rho/(1.0 + 2*rho*epsize[i])*exp(-2*rho*tk[i]);
    tmp -= 2.0*exp(-(tk[i+1] - tk[i])/epsize[i] - 2.0*rho*tk[i] );
    tmp -= 2.0*rho*epsize[i]/(1.0 - epsize[i]*2.0*rho)*exp(-2.0*rho*tk[i]-2.0*(tk[i+1]-tk[i])/epsize[i]);
    tmp += 2.0/(1.0-epsize[i]*2.0*rho)/(1.0 + 2.0*rho)*exp(-rho*tk[i+1]-(tk[i+1]-tk[i])/epsize[i]);
    P4[i] *= tmp;
  }
  unsigned i = tk_l - 1;
  P4[i] = 2.0*rho/(1.0 + 2.0*rho*epsize[i])*exp(-2.0*rho*tk[i]);
}
	
  
void ComputeP7(double *tk,int tk_l,double *P7,double *epsize,double rho){
  for (unsigned i = 0; i < tk_l - 1; i++){
    P7[i] = 1.0 - exp(-(tk[i+1]-tk[i])*2.0*rho) - exp(-tk[i]*2.0*rho);
    P7[i] -= epsize[i]*2*rho/(1 - epsize[i]*2.0*rho)*exp(-(tk[i+1]-tk[i])/epsize[i]-tk[i]*2.0*rho);
    P7[i] += 1.0/(1.0 - epsize[i]*2.0*rho)*exp(-tk[i]*2.0*rho);
  }
  unsigned i = tk_l - 1;
  P7[i] = 1.0 - exp(-2.0*rho*tk[i]);
}
  
void ComputeP0(int tk_l,double *P0,double *P5){ //probability P(T > i)
  P0[0] = P5[0];
  for (unsigned i = 1; i < tk_l; i++)
    P0[i] = P0[i-1]*P5[i];
}
