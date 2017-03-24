
void ComputeP11(unsigned numWind,int tk_l,double **P,double **PP,double **fw,double **bw,double *stationary){
  for (unsigned i = 0; i < tk_l; i++){
    PP[1][i] = 0;
    for (unsigned l = 1; l < numWind; l++)
      PP[1][i] += fw[i][l]*P[1][i]*bw[l+1][i]/stationary[i];//NOTE: In appendix of K.H. paper it seems to be an extra emission probability for site l+1, it is already inside bw[]
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

void ComputeP33(unsigned numWind,int tk_l,double **P,double **PP,double **fw,double **bw,double *stationary){
  double R1[tk_l];
  for (unsigned i = 0; i < tk_l; i++)
    PP[3][i] = 0;
  for (unsigned l = 1; l < numWind; l++){
    R1[tk_l - 1] = 0;
    for (int i = tk_l - 2; i >= 0 ; i--)
      R1[i] = R1[i+1] + fw[i+1][l];
    for (unsigned i = 0; i < tk_l - 1; i++)
      PP[3][i] += R1[i]*P[3][i]*bw[i][l]/stationary[i];
  }
}

void ComputeP44(unsigned numWind,int tk_l,double **P,double **PP,double **fw,double **bw,double *stationary){
  for (unsigned i = 0; i < tk_l; i++){
    PP[4][i] = 0;
    for (unsigned l = 1; l < numWind; l++)
      PP[4][i] += fw[i][l]*P[4][i]*bw[i][l+1]/stationary[i];
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

