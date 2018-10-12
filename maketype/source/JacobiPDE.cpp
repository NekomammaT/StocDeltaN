#include "JacobiPDE.hpp"

// ---------------------- sample ------------------------

double JacobiPDE::H(vector<double> &X, vector<double> &P)
{
  double rho = V(X);

  for (int I=0; I<X.size(); I++) {
    for (int J=0; J<X.size(); J++) {
      rho += 1./2*inversemetric(X,I,J)*P[I]*P[J];
    }
  }

  return sqrt(rho/3.);
}

double JacobiPDE::V(vector<double> &X)
{
  double m1 = 0.01;
  double m2 = 0.1;
  
  return 1./2*m1*m1*X[0]*X[0] + 1./2*m2*m2*X[1]*X[1];
}

double JacobiPDE::VI(vector<double> &X, int I)
{
  double m1 = 0.01;
  double m2 = 0.1;

  if (I == 0) {
    return m1*m1*X[0];
  } else {
    return m2*m2*X[1];
  }
}

double JacobiPDE::metric(vector<double> &X, int I, int J)
{
  double MM = 1e-3;

  if (I == 0 && J == 0) {
    return 1 + 2*X[1]*X[1]/MM/MM;
  } else if (I == 1 && J == 1) {
    return 1;
  } else {
    return 0;
  }
}

double JacobiPDE::inversemetric(vector<double> &X, int I, int J)
{
  double MM = 1e-3;

  if (I == 0 && J == 0) {
    return 1./(1+2*X[1]*X[1]/MM/MM);
  } else if (I == 1 && J == 1) {
    return 1;
  } else {
    return 0;
  }
}

double JacobiPDE::affine(vector<double> &X, int I, int J, int K)
{
  double MM = 1e-3;

  if (I == 0 && ((J == 0 && K == 1) || (J == 1 && K == 0))) {
    return 2*X[1] / (2*X[1]*X[1] + MM*MM);
  } else if (I == 1 && J == 0 && K == 0) {
    return -2*X[1]/MM/MM;
  } else {
    return 0;
  }
}

double JacobiPDE::Dphi(vector<double> &X, vector<double> &P, int I)
{
  double Dphi = 0;

  for (int J=0; J<dim; J++) {
    Dphi += inversemetric(X,I,J)*P[J]/H(X,P);
  }

  return Dphi;
}

double JacobiPDE::Dpi(vector<double> &X, vector<double> &P, int I)
{
  double Hubble = H(X,P);
  double Dpi = -3*P[I]-VI(X,I)/Hubble;

  for (int J=0; J<dim; J++) {
    for (int K=0; K<dim; K++) {
      for (int L=0; L<dim; L++) {
	Dpi += affine(X,K,I,J)*inversemetric(X,J,L)*P[K]*P[L]/Hubble;
      }
    }
  }

  return Dpi;
}

double JacobiPDE::Dphiphi(vector<double> &X, vector<double> &P, int I, int J)
{
  return H(X,P)*H(X,P)/4./M_PI/M_PI * inversemetric(X,I,J);
}

double JacobiPDE::Dphipi(vector<double> &X, vector<double> &P, int I, int J)
{
  double Dphipi = 0;

  for (int K=0; K<dim; K++) {
    for (int L=0; L<dim; L++) {
      Dphipi += affine(X,K,J,L)*P[K]*Dphiphi(X,P,I,L);
    }
  }

  return Dphipi;
}

double JacobiPDE::Dpipi(vector<double> &X, vector<double> &P, int I, int J)
{
  double Dpipi = 0;

  for (int K=0; K<dim; K++) {
    for (int L=0; L<dim; L++) {
      for (int M=0; M<dim; M++) {
	for (int N=0; N<dim; N++) {
	  Dpipi += affine(X,K,I,L)*P[K]*affine(X,M,J,N)*P[M]*Dphiphi(X,P,L,N);
	}
      }
    }
  }

  return Dpipi;
}

// -------------------------------------------------------



JacobiPDE::JacobiPDE(vector< vector<double> > Site[], double Rhoc)
{
  srand((unsigned)time(NULL));

#ifdef _OPENMP
  cout<<"OpenMP : Enabled (Max # of threads = "<<omp_get_max_threads()<<")"<<endl;
#endif

  rhoc = Rhoc;
  site[0] = Site[0];
  site[1] = Site[1];
  hI[0] = site[0];
  hI[1] = site[1];
  dim = site[0].size();
  siteNo[0] = vector<int>(dim,0);
  siteNo[1] = vector<int>(dim,0);
  ind[0] = vector<int>(dim,0);
  ind[1] = vector<int>(dim,0);
  FPoint[0] = vector<double>(dim,0);
  FPoint[1] = vector<double>(dim,0);

  for (int xp=0; xp<2; xp++) {
    for (int I=0; I<dim; I++) {
      siteNo[xp][I] = site[xp][I].size();
      
      for (int s=0; s<site[xp][I].size()-1; s++) {
	hI[xp][I][s] = site[xp][I][s+1] - site[xp][I][s];
      }
    }
  }

  volume = 1, pvol = 1;
  for (int I=0; I<dim; I++) {
    volume *= siteNo[0][I]*siteNo[1][I];
    pvol *= siteNo[1][I];
  }

  f1 = vector<double>(volume,0);
  f1_next = f1;
  g2 = f1;
  g2_next = f1;
  Omega = vector<bool>(volume,true);
  
  for (int number=0; number<volume; number++) {
    for (int I=0; I<dim; I++) {
      FPoint[0][I] = No2X(number,I);
      FPoint[1][I] = No2P(number,I);
    }
    
    if (V(FPoint[0]) < rhoc) {
      Omega[number] = false;
      f1[number] = 0;
      g2[number] = 0;
    } else {
      f1[number] = rand()%10;
      g2[number] = (rand()%10)/10.;
    }
  }
  
  // ---------------- reflecting b.c. ----------------

  numXp = vector< vector<int> >(volume, vector<int>(dim,0));
  numXm = numXp;
  numPp = numXp;
  numPm = numXp;

  numXXpp
    = vector< vector< vector<int> > >(volume, vector< vector<int> >(dim, vector<int>(dim,0)));
  numXXpm = numXXpp;
  numXXmm = numXXpp;
  numXPpp = numXXpp;
  numXPpm = numXXpp;
  numXPmp = numXXpp;
  numXPmm = numXXpp;
  numPPpp = numXXpp;
  numPPpm = numXXpp;
  numPPmm = numXXpp;

  hp[0] = vector< vector<double> >(volume, vector<double>(dim,0));
  hp[1] = hp[0];
  hm[0] = hp[0];
  hm[1] = hp[0];

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int number=0; number<volume; number++) {
    vector<int> ind0[2],indXp[2],indXm[2],indPp[2],indPm[2],
      indXXpp[2],indXXpm[2],indXXmm[2],indXPpp[2],indXPpm[2],indXPmp[2],indXPmm[2],
      indPPpp[2],indPPpm[2],indPPmm[2];

    ind0[0] = ind[0];
    ind0[1] = ind[1];

    for (int I=0; I<dim; I++) {
      ind0[0][I] = No2XInd(number,I);
      ind0[1][I] = No2PInd(number,I);
    }

    for (int I=0; I<dim; I++) {
      for (int xp=0; xp<=1; xp++) {
	indXp[xp] = ind0[xp];
	indXm[xp] = ind0[xp];
	indPp[xp] = ind0[xp];
	indPm[xp] = ind0[xp];
      }

      for (int J=0; J<dim; J++) {
	if (J != I) {
	  for (int xp=0; xp<=1; xp++) {
	    indXXpp[xp] = ind0[xp];
	    indXXpm[xp] = ind0[xp];
	    indXXmm[xp] = ind0[xp];
	    indXPpp[xp] = ind0[xp];
	    indXPpm[xp] = ind0[xp];
	    indXPmp[xp] = ind0[xp];
	    indXPmm[xp] = ind0[xp];
	    indPPpp[xp] = ind0[xp];
	    indPPpm[xp] = ind0[xp];
	    indPPmm[xp] = ind0[xp];
	  }
	  if (ind0[0][I] == 0) {
	    indXXmm[0][I]++;
	    indXPmp[0][I]++;
	    indXPmm[0][I]++;
	  } else {
	    indXXmm[0][I]--;
	    indXPmp[0][I]--;
	    indXPmm[0][I]--;
	  }
	  if (ind0[0][I] == siteNo[0][I]-1) {
	    indXXpp[0][I]--;
	    indXXpm[0][I]--;
	    indXPpp[0][I]--;
	    indXPpm[0][I]--;
	  } else {
	    indXXpp[0][I]++;
	    indXXpm[0][I]++;
	    indXPpp[0][I]++;
	    indXPpm[0][I]++;
	  }
	  if (ind0[0][J] == 0) {
	    indXXpm[0][J]++;
	    indXXmm[0][J]++;
	  } else {
	    indXXpm[0][J]--;
	    indXXmm[0][J]--;
	  }
	  if (ind0[0][J] == siteNo[0][J]-1) {
	    indXXpp[0][J]--;
	  } else {
	    indXXpp[0][J]++;
	  }
	  if (ind0[1][I] == 0) {
	    indPPmm[1][I]++;
	  } else {
	    indPPmm[1][I]--;
	  }
	  if (ind0[1][I] == siteNo[1][I]-1) {
	    indPPpp[1][I]--;
	    indPPpm[1][I]--;
	  } else {
	    indPPpp[1][I]++;
	    indPPpm[1][I]++;
	  }
	  if (ind0[1][J] == 0) {
	    indXPpm[1][J]++;
	    indXPmm[1][J]++;
	    indPPpm[1][J]++;
	    indPPmm[1][J]++;
	  } else {
	    indXPpm[1][J]--;
	    indXPmm[1][J]--;
	    indPPpm[1][J]--;
	    indPPmm[1][J]--;
	  }
	  if (ind0[1][J] == siteNo[1][J]-1) {
	    indXPpp[1][J]--;
	    indXPmp[1][J]--;
	    indPPpp[1][J]--;
	  } else {
	    indXPpp[1][J]++;
	    indXPmp[1][J]++;
	    indPPpp[1][J]++;
	  }

	  numXXpp[number][I][J] = Ind2No(indXXpp);
	  numXXpm[number][I][J] = Ind2No(indXXpm);
	  numXXmm[number][I][J] = Ind2No(indXXmm);
	  numXPpp[number][I][J] = Ind2No(indXPpp);
	  numXPpm[number][I][J] = Ind2No(indXPpm);
	  numXPmp[number][I][J] = Ind2No(indXPmp);
	  numXPmm[number][I][J] = Ind2No(indXPmm);
	  numPPpp[number][I][J] = Ind2No(indPPpp);
	  numPPpm[number][I][J] = Ind2No(indPPpm);
	  numPPmm[number][I][J] = Ind2No(indPPmm);
	}
      }

      if (ind0[0][I] == 0) {
	indXm[0][I]++;
	hm[0][number][I] = hI[0][I][ind0[0][I]];
      } else {
	indXm[0][I]--;
	hm[0][number][I] = hI[0][I][ind0[0][I]-1];
      }
      if (ind0[0][I] == siteNo[0][I]-1) {
	indXp[0][I]--;
	hp[0][number][I] = hI[0][I][ind0[0][I]-1];
      } else {
	indXp[0][I]++;
	hp[0][number][I] = hI[0][I][ind0[0][I]];
      }
      if (ind0[1][I] == 0) {
	indPm[1][I]++;
	hm[1][number][I] = hI[1][I][ind0[1][I]];
      } else {
	indPm[0][I]--;
	hm[1][number][I] = hI[1][I][ind0[1][I]-1];
      }
      if (ind0[1][I] == siteNo[1][I]-1) {
	indPp[1][I]--;
	hp[1][number][I] = hI[1][I][ind0[1][I]-1];
      } else {
	indPp[1][I]++;
	hp[1][number][I] = hI[1][I][ind0[1][I]];
      }

      numXp[number][I] = Ind2No(indXp);
      numXm[number][I] = Ind2No(indXm);
      numPp[number][I] = Ind2No(indPp);
      numPm[number][I] = Ind2No(indPm);
    }
  }
}

double JacobiPDE::PDE(int number, int n)
{
  vector<double> FPoint0[2];

  FPoint0[0] = FPoint[0];
  FPoint0[1] = FPoint[1];
  
  for (int I=0; I<dim; I++) {
    FPoint0[0][I] = No2X(number,I);
    FPoint0[1][I] = No2P(number,I);
  }
  
  double uXIm, uXIp, uPIm, uPIp, uXXmm, uXXmp, uXXpm, uXXpp,
    uXPmm, uXPmp, uXPpm, uXPpp, uPPmm, uPPmp, uPPpm, uPPpp,
    hXIp, hXIm, hPIp, hPIm, hXJp, hXJm, hPJp, hPJm;
  double uij = 0, coeff = 0;

  if (n == 1) {
    uij -= 1;
  } // Cn
  
  for (int I=0; I<dim; I++) {
    hXIp = hp[0][number][I];
    hXIm = hm[0][number][I];
    hPIp = hp[1][number][I];
    hPIm = hm[1][number][I];

    if (n == 1) {
      uXIm = f1[numXm[number][I]];
      uXIp = f1[numXp[number][I]];
      uPIm = f1[numPm[number][I]];
      uPIp = f1[numPp[number][I]];
    } else if (n == -2) {
      uXIm = g2[numXm[number][I]];
      uXIp = g2[numXp[number][I]];
      uPIm = g2[numPm[number][I]];
      uPIp = g2[numPp[number][I]];
    }

    if (Dphi(FPoint0[0],FPoint0[1],I) < 0 ) {
      uij += Dphi(FPoint0[0],FPoint0[1],I)/hXIm * uXIm;
      coeff += Dphi(FPoint0[0],FPoint0[1],I)/hXIm;
    } else {
      uij -= Dphi(FPoint0[0],FPoint0[1],I)/hXIp * uXIp;
      coeff -= Dphi(FPoint0[0],FPoint0[1],I)/hXIp;
    }

    if (Dpi(FPoint0[0],FPoint0[1],I) < 0) {
      uij += Dpi(FPoint0[0],FPoint0[1],I)/hPIm * uPIm;
      coeff += Dpi(FPoint0[0],FPoint0[1],I)/hPIm;
    } else {
      uij -= Dpi(FPoint0[0],FPoint0[1],I)/hPIp * uPIp;
      coeff -= Dpi(FPoint0[0],FPoint0[1],I)/hPIp;
    }

    uij -= Dphiphi(FPoint0[0],FPoint0[1],I,I)
      *(uXIp*hXIm+uXIm*hXIp)/hXIp/hXIm/(hXIm+hXIp)
      + Dpipi(FPoint0[0],FPoint0[1],I,I)
      *(uPIp*hPIm+uPIm*hPIp)/hPIp/hPIm/(hPIm+hPIp);
    coeff -= Dphiphi(FPoint0[0],FPoint0[1],I,I)/hXIp/hXIm
      + Dpipi(FPoint0[0],FPoint0[1],I,I)/hPIp/hPIm;
    
    
    for (int J=0; J<dim; J++) {
      hXJp = hp[0][number][J];
      hXJm = hm[0][number][J];
      hPJp = hp[1][number][J];
      hPJm = hm[1][number][J];

      if (J != I) {
	if (n == 1) {
	  uXXmm = f1[numXXmm[number][I][J]];
	  uXXmp = f1[numXXpm[number][J][I]];
	  uXXpm = f1[numXXpm[number][I][J]];
	  uXXpp = f1[numXXpp[number][I][J]];
	  
	  uXPmm = f1[numXPmm[number][I][J]];
	  uXPmp = f1[numXPmp[number][I][J]];
	  uXPpm = f1[numXPpm[number][I][J]];
	  uXPpp = f1[numXPpp[number][I][J]];
	  
	  uPPmm = f1[numPPmm[number][I][J]];
	  uPPmp = f1[numPPpm[number][J][I]];
	  uPPpm = f1[numPPpm[number][I][J]];
	  uPPpp = f1[numPPpp[number][I][J]];
	} else if (n == -2) {
	  uXXmm = g2[numXXmm[number][I][J]];
	  uXXmp = g2[numXXpm[number][J][I]];
	  uXXpm = g2[numXXpm[number][I][J]];
	  uXXpp = g2[numXXpp[number][I][J]];
	  
	  uXPmm = g2[numXPmm[number][I][J]];
	  uXPmp = g2[numXPmp[number][I][J]];
	  uXPpm = g2[numXPpm[number][I][J]];
	  uXPpp = g2[numXPpp[number][I][J]];
	  
	  uPPmm = g2[numPPmm[number][I][J]];
	  uPPmp = g2[numPPpm[number][J][I]];
	  uPPpm = g2[numPPpm[number][I][J]];
	  uPPpp = g2[numPPpp[number][I][J]];
	}
	
	uij -=
	  1./2*Dphiphi(FPoint0[0],FPoint0[1],I,J)*(uXXpp-uXXpm-uXXmp+uXXmm)
	  /(hXIm+hXIp)/(hXJm+hXJp)
	  + Dphipi(FPoint0[0],FPoint0[1],I,J)*(uXPpp-uXPpm-uXPmp+uXPmm)
	  /(hXIm+hXIp)/(hPJm+hPJp)
	  + 1./2*Dpipi(FPoint0[0],FPoint0[1],I,J)*(uPPpp-uPPpm-uPPmp+uPPmm)
	  /(hPIm+hPIp)/(hPJm+hPJp);
      }
      
      if (n == -2) {
	uij -=
	  Dphiphi(FPoint0[0],FPoint0[1],I,J)*(f1[numXp[number][I]]-f1[numXm[number][I]])
	  *(f1[numXp[number][J]]-f1[numXm[number][J]])/(hXIm+hXIp)/(hXJm+hXJp)
	  + 2*Dphipi(FPoint0[0],FPoint0[1],I,J)*(f1[numXp[number][I]]-f1[numXm[number][I]])
	  *(f1[numPp[number][J]]-f1[numPm[number][J]])/(hXIm+hXIp)/(hPJm+hPJp)
	  + Dpipi(FPoint0[0],FPoint0[1],I,J)*(f1[numPp[number][I]]-f1[numPm[number][I]])
	  *(f1[numPp[number][J]]-f1[numPm[number][J]])/(hPIm+hPIp)/(hPJm+hPJp);
      } // Cn
    }
  }

  uij /= coeff;

  return uij;
}

void JacobiPDE::PDE_solve(int maxstep, double tol, int n)
{
  double unext, u_norm, err;
  
  for (int step=0; step<maxstep; step++) {
    u_norm = 0;
    err = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:u_norm, err)
#endif
    {
#ifdef _OPENMP
#pragma omp for
#endif
      for (int i=0; i<volume; i++) {
	if (Omega[i]) {
	  if (n == 1) {
	    u_norm += f1[i]*f1[i];
	    f1_next[i] = PDE(i,n);
	    err += (f1[i]-f1_next[i])*(f1[i]-f1_next[i]);
	  } else if (n == -2) {
	    u_norm += g2[i]*g2[i];
	    g2_next[i] = PDE(i,n);
	    err += (g2[i]-g2_next[i])*(g2[i]-g2_next[i]);
	  }
	}
      }
    }
    
    if (n == 1) {
      f1 = f1_next;
    } else if (n == -2) {
      g2 = g2_next;
    }
    
    err = sqrt(err)/sqrt(u_norm);

    if (n == 1) {
      cout << "\rerr for M1: " << setw(11) << left << err << "  step: " << step << flush;
    } else if (n == -2) {
      cout << "\rerr for C2: " << setw(11) << left << err << "  step: " << step << flush;
    }
    
    if (err < tol) {
      break;
    }
  }

  cout << endl;
}

int JacobiPDE::Ind2No(vector<int> index[])
{
  int number = 0;
  int xtemp, ptemp;

  for (int I=0; I<dim; I++) {
    xtemp = index[0][dim-1-I] * pvol;
    ptemp = index[1][dim-1-I];
    for (int J=0; J<I; J++) {
      xtemp *= siteNo[0][dim-1-J];
      ptemp *= siteNo[1][dim-1-J];
    }
    number += xtemp + ptemp;
  }
  
  return number;
}

int JacobiPDE::No2XInd(int number, int I)
{
  int index = number / pvol;

  for (int J=dim-1; J>I; J--) {
    index /= siteNo[0][J];
  }

  if (I > 0) {
    index %= siteNo[0][I];
  }

  return index;
}

int JacobiPDE::No2PInd(int number, int I)
{
  int index = number;

  for (int J=dim-1; J>I; J--) {
    index /= siteNo[1][J];
  }

  return index % siteNo[1][I];
}

double JacobiPDE::No2X(int number, int I)
{
  return site[0][I][No2XInd(number,I)];
}

double JacobiPDE::No2P(int number, int I)
{
  return site[1][I][No2PInd(number,I)];
}

int JacobiPDE::ceilX(vector<double> &X, int I)
{
  int xind = 0;

  while (site[0][I][xind] <= X[I]) {
    xind++;
  }

  return xind;
}

int JacobiPDE::ceilP(vector<double> &P, int I)
{
  int pind = 0;

  while (site[1][I][pind] <= P[I]) {
    pind++;
  }

  return pind;
}

double JacobiPDE::Interpolation_f(vector<double> &X, vector<double> &P, vector<double> &f) 
{
  int pmcheckX, pmcheckP;
  double weight, intf = 0;
  vector< vector<int> > near_site[2];
  vector<double> fdata;
  
  for (int xnum=0; xnum<pow(2,dim); xnum++) {
    for (int pnum=0; pnum<pow(2,dim); pnum++) {
      weight = 1;
      
      for (int I=0; I<dim; I++) {
	pmcheckX = xnum;
	pmcheckP = pnum;
      
	for (int J=0; J<I; J++) {
	  pmcheckX /= 2;
	  pmcheckP /= 2;
	}
	
	if (I < dim-1) {
	  pmcheckX %= 2;
	  pmcheckP %= 2;
	}

	if (pmcheckX == 1) {
	  ind[0][I] = ceilX(X,I);
	  weight *= (X[I]-site[0][I][ceilX(X,I)-1])/hI[0][I][ceilX(X,I)-1];
	} else {
	  ind[0][I] = ceilX(X,I)-1;
	  weight *= (site[0][I][ceilX(X,I)]-X[I])/hI[0][I][ceilX(X,I)-1];
	}
	if (pmcheckP == 1) {
	  ind[1][I] = ceilP(P,I);
	  weight *= (P[I]-site[1][I][ceilP(P,I)-1])/hI[1][I][ceilP(P,I)-1];
	} else {
	  ind[1][I] = ceilP(P,I)-1;
	  weight *= (site[1][I][ceilP(P,I)]-P[I])/hI[1][I][ceilP(P,I)-1];
	}
      }

      intf += f[Ind2No(ind)]*weight;
    }
  }

  return intf;
}

double JacobiPDE::return_f1(vector<int> index[])
{
  return f1[Ind2No(index)];
}

double JacobiPDE::return_g2(vector<int> index[])
{
  return g2[Ind2No(index)];
}

void JacobiPDE::export_fg(string filename)
{
  ofstream ofs(filename);

  for (int number=0; number<volume; number++) {
    for (int I=0; I<dim; I++) {
      ofs << No2X(number,I) << ' ';
    }
    for (int I=0; I<dim; I++) {
      ofs << No2P(number,I) << ' ';
    }
    ofs << f1[number] << ' '
	<< g2[number] << endl;
  }
}
