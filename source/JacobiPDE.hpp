#ifndef INCLUDED_JacobiPDE_hpp_
#define INCLUDED_JacobiPDE_hpp_

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <omp.h>

using namespace std;

class JacobiPDE
{
protected:
  vector<double> f1, f1_next, g2, g2_next; // argument: number
  vector< vector<double> > site[2], hI[2]; // argument: x or p, I, siteNo
  vector<double> FPoint[2], exportdx[2]; // argument: x or p, I
  vector<bool> Omega; // argument: number
  vector<int> ind[2]; // argument: x or p, I
  vector<int> siteNo[2]; // argument: x or p, I
  double rhoc, err;
  int dim, volume, pvol;
  
public:
  JacobiPDE(){}
  JacobiPDE(vector< vector<double> > Site[], double Rhoc, vector<double> Exportdx[]);
  void init_fn();
  double PDE(int number, int n);
  void PDE_solve(int maxstep, double tol, int n);
  int Ind2No(vector<int> index[]);
  int No2XInd(int number, int I);
  int No2PInd(int number, int I);
  double No2X(int number, int I);
  double No2P(int number, int I);
  int ceilX(vector<double> &X, int I);
  int ceilP(vector<double> &X, int I);
  double Interpolation_f(vector<double> &X, vector<double> &P, vector<double> &f);
  double return_f1(vector<int> index[]);
  double return_g2(vector<int> index[]);
  double return_err();
  void export_fg(string filename);
  double H(vector<double> &X, vector<double> &P);
  virtual double V(vector<double> &X);
  virtual double VI(vector<double> &X, int I);
  virtual double metric(vector<double> &X, int I, int J);
  virtual double inversemetric(vector<double> &X, int I, int J);
  virtual double affine(vector<double> &X, int I, int J, int K);
  virtual double Dphi(vector<double> &X, vector<double> &P, int I);
  virtual double Dpi(vector<double> &X, vector<double> &P, int I);
  virtual double Dphiphi(vector<double> &X, vector<double> &P, int I, int J);
  virtual double Dphipi(vector<double> &X, vector<double> &P, int I, int J);
  virtual double Dpipi(vector<double> &X, vector<double> &P, int I, int J);
};



// ---------------------- sample ------------------------

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



JacobiPDE::JacobiPDE(vector< vector<double> > Site[], double Rhoc, vector<double> Exportdx[])
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
  exportdx[0] = Exportdx[0];
  exportdx[1] = Exportdx[1];

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
}

void JacobiPDE::init_fn()
{
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
}

double JacobiPDE::PDE(int number, int n)
{
  vector<int> ind0[2];
  vector<double> FPoint0[2];

  ind0[0] = ind[0];
  ind0[1] = ind[1];
  FPoint0[0] = FPoint[0];
  FPoint0[1] = FPoint[1];
  
  for (int I=0; I<dim; I++) {
    ind0[0][I] = No2XInd(number,I);
    ind0[1][I] = No2PInd(number,I);
    FPoint0[0][I] = No2X(number,I);
    FPoint0[1][I] = No2P(number,I);
  }
  
  vector<int> indXIm[2], indXIp[2], indXJm[2], indXJp[2],
    indPIm[2], indPIp[2], indPJm[2], indPJp[2],
    indXXmm[2], indXXmp[2], indXXpm[2], indXXpp[2],
    indXPmm[2], indXPmp[2], indXPpm[2], indXPpp[2],
    indPPmm[2], indPPmp[2], indPPpm[2], indPPpp[2];
  double uXIm, uXIp, uPIm, uPIp, uXXmm, uXXmp, uXXpm, uXXpp,
    uXPmm, uXPmp, uXPpm, uXPpp, uPPmm, uPPmp, uPPpm, uPPpp,
    hXIp, hXIm, hPIp, hPIm, hXJp, hXJm, hPJp, hPJm;
  double uij = 0, coeff = 0;

  if (n == 1) {
    uij -= 1;
  } // Cn
  
  for (int I=0; I<dim; I++) {
    for (int xp=0; xp<=1; xp++) {
      indXIm[xp] = ind0[xp];
      indXIp[xp] = ind0[xp];
      indPIm[xp] = ind0[xp];
      indPIp[xp] = ind0[xp];
    }
    
    // -------------- reflecting b.c. -----------------
    
    if (indXIm[0][I] == 0) {
      indXIm[0][I]++;
      hXIm = hI[0][I][ind0[0][I]];
    } else {
      indXIm[0][I]--;
      hXIm = hI[0][I][ind0[0][I]-1];
    }
    if (indPIm[1][I] == 0) {
      indPIm[1][I]++;
      hPIm = hI[1][I][ind0[1][I]];
    } else {
      indPIm[1][I]--;
      hPIm = hI[1][I][ind0[1][I]-1];
    }

    if (indXIp[0][I] == siteNo[0][I]-1) {
      indXIp[0][I]--;
      hXIp = hI[0][I][ind0[0][I]-1];
    } else {
      indXIp[0][I]++;
      hXIp = hI[0][I][ind0[0][I]];
    }
    if (indPIp[1][I] == siteNo[1][I]-1) {
      indPIp[1][I]--;
      hPIp = hI[1][I][ind0[1][I]-1];
    } else {
      indPIp[1][I]++;
      hPIp = hI[1][I][ind0[1][I]];
    }

    // ----------------------------------------------

    if (n == 1) {
      uXIm = f1[Ind2No(indXIm)];
      uXIp = f1[Ind2No(indXIp)];
      uPIm = f1[Ind2No(indPIm)];
      uPIp = f1[Ind2No(indPIp)];
    } else if (n == -2) {
      uXIm = g2[Ind2No(indXIm)];
      uXIp = g2[Ind2No(indXIp)];
      uPIm = g2[Ind2No(indPIm)];
      uPIp = g2[Ind2No(indPIp)];
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
      for (int xp=0; xp<=1; xp++) {
	indXJm[xp] = ind0[xp];
	indXJp[xp] = ind0[xp];
	indPJm[xp] = ind0[xp];
	indPJp[xp] = ind0[xp];
	indXXmm[xp] = ind0[xp];
	indXXmp[xp] = ind0[xp];
	indXXpm[xp] = ind0[xp];
	indXXpp[xp] = ind0[xp];
	indXPmm[xp] = ind0[xp];
	indXPmp[xp] = ind0[xp];
	indXPpm[xp] = ind0[xp];
	indXPpp[xp] = ind0[xp];
	indPPmm[xp] = ind0[xp];
	indPPmp[xp] = ind0[xp];
	indPPpm[xp] = ind0[xp];
	indPPpp[xp] = ind0[xp];
      }
	
      // -------------------- reflecting b.c. --------------------
      
      if (indXJm[0][J] == 0) {
	indXJm[0][J]++;
	hXJm = hI[0][J][ind0[0][J]];
      } else {
	indXJm[0][J]--;
	hXJm = hI[0][J][ind0[0][J]-1];
      }
      if (indPJm[1][J] == 0) {
	indPJm[1][J]++;
	hPJm = hI[1][J][ind0[1][J]];
      } else {
	indPJm[1][J]--;
	hPJm = hI[1][J][ind0[1][J]-1];
      }
      if (indXJp[0][J] == siteNo[0][J]-1) {
	indXJp[0][J]--;
	hXJp = hI[0][J][ind0[0][J]-1];
      } else {
	indXJp[0][J]++;
	hXJp = hI[0][J][ind0[0][J]];
      }
      if (indPJp[1][J] == siteNo[1][J]-1) {
	indPJp[1][J]--;
	hPJp = hI[1][J][ind0[1][J]-1];
      } else {
	indPJp[1][J]++;
	hPJp = hI[1][J][ind0[1][J]];
      }

      
      if (indXXmm[0][I] == 0) {
	indXXmm[0][I]++;
      } else {
	indXXmm[0][I]--;
      }
      if (indXXmp[0][I] == 0) {
	indXXmp[0][I]++;
      } else {
	indXXmp[0][I]--;
      }
      if (indXXpm[0][I] == siteNo[0][I]-1) {
	indXXpm[0][I]--;
      } else {
	indXXpm[0][I]++;
      }
      if (indXXpp[0][I] == siteNo[0][I]-1) {
	indXXpp[0][I]--;
      } else {
	indXXpp[0][I]++;
      }
      
      if (indXXmm[0][J] == 0) {
	indXXmm[0][J]++;
      } else {
	indXXmm[0][J]--;
      }
      if (indXXpm[0][J] == 0) {
	indXXpm[0][J]++;
      } else {
	indXXpm[0][J]--;
      }
      if (indXXmp[0][J] == siteNo[0][J]-1) {
	indXXmp[0][J]--;
      } else {
	indXXmp[0][J]++;
      }
      if (indXXpp[0][J] == siteNo[0][J]-1) {
	indXXpp[0][J]--;
      } else {
	indXXpp[0][J]++;
      }


      if (indXPmm[0][I] == 0) {
	indXPmm[0][I]++;
      } else {
	indXPmm[0][I]--;
      }
      if (indXPmp[0][I] == 0) {
	indXPmp[0][I]++;
      } else {
	indXPmp[0][I]--;
      }
      if (indXPpm[0][I] == siteNo[0][I]-1) {
	indXPpm[0][I]--;
      } else {
	indXPpm[0][I]++;
      }
      if (indXPpp[0][I] == siteNo[0][I]-1) {
	indXPpp[0][I]--;
      } else {
	indXPpp[0][I]++;
      }
      
      if (indXPmm[1][J] == 0) {
	indXPmm[1][J]++;
      } else {
	indXPmm[1][J]--;
      }
      if (indXPpm[1][J] == 0) {
	indXPpm[1][J]++;
      } else {
	indXPpm[1][J]--;
      }
      if (indXPmp[1][J] == siteNo[1][J]-1) {
	indXPmp[1][J]--;
      } else {
	indXPmp[1][J]++;
      }
      if (indXPpp[1][J] == siteNo[1][J]-1) {
	indXPpp[1][J]--;
      } else {
	indXPpp[1][J]++;
      }


      if (indPPmm[1][I] == 0) {
	indPPmm[1][I]++;
      } else {
	indPPmm[1][I]--;
      }
      if (indPPmp[1][I] == 0) {
	indPPmp[1][I]++;
      } else {
	indPPmp[1][I]--;
      }
      if (indPPpm[1][I] == siteNo[1][I]-1) {
	indPPpm[1][I]--;
      } else {
	indPPpm[1][I]++;
      }
      if (indPPpp[1][I] == siteNo[1][I]-1) {
	indPPpp[1][I]--;
      } else {
	indPPpp[1][I]++;
      }
      
      if (indPPmm[1][J] == 0) {
	indPPmm[1][J]++;
      } else {
	indPPmm[1][J]--;
      }
      if (indPPpm[1][J] == 0) {
	indPPpm[1][J]++;
      } else {
	indPPpm[1][J]--;
      }
      if (indPPmp[1][J] == siteNo[1][J]-1) {
	indPPmp[1][J]--;
      } else {
	indPPmp[1][J]++;
      }
      if (indPPpp[1][J] == siteNo[1][J]-1) {
	indPPpp[1][J]--;
      } else {
	indPPpp[1][J]++;
      }
      
      // ----------------------------------

      if (n == 1) {
	uXXmm = f1[Ind2No(indXXmm)];
	uXXmp = f1[Ind2No(indXXmp)];
	uXXpm = f1[Ind2No(indXXpm)];
	uXXpp = f1[Ind2No(indXXpp)];

	uXPmm = f1[Ind2No(indXPmm)];
	uXPmp = f1[Ind2No(indXPmp)];
	uXPpm = f1[Ind2No(indXPpm)];
	uXPpp = f1[Ind2No(indXPpp)];

	uPPmm = f1[Ind2No(indPPmm)];
	uPPmp = f1[Ind2No(indPPmp)];
	uPPpm = f1[Ind2No(indPPpm)];
	uPPpp = f1[Ind2No(indPPpp)];
      } else if (n == -2) {
	uXXmm = g2[Ind2No(indXXmm)];
	uXXmp = g2[Ind2No(indXXmp)];
	uXXpm = g2[Ind2No(indXXpm)];
	uXXpp = g2[Ind2No(indXXpp)];

	uXPmm = g2[Ind2No(indXPmm)];
	uXPmp = g2[Ind2No(indXPmp)];
	uXPpm = g2[Ind2No(indXPpm)];
	uXPpp = g2[Ind2No(indXPpp)];

	uPPmm = g2[Ind2No(indPPmm)];
	uPPmp = g2[Ind2No(indPPmp)];
	uPPpm = g2[Ind2No(indPPpm)];
	uPPpp = g2[Ind2No(indPPpp)];
      }
             
      if (J != I) {
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
	  Dphiphi(FPoint0[0],FPoint0[1],I,J)*(f1[Ind2No(indXIp)]-f1[Ind2No(indXIm)])
	  *(f1[Ind2No(indXJp)]-f1[Ind2No(indXJm)])/(hXIm+hXIp)/(hXJm+hXJp)
	  + 2*Dphipi(FPoint0[0],FPoint0[1],I,J)*(f1[Ind2No(indXIp)]-f1[Ind2No(indXIm)])
	  *(f1[Ind2No(indPJp)]-f1[Ind2No(indPJm)])/(hXIm+hXIp)/(hPJm+hPJp)
	  + Dpipi(FPoint0[0],FPoint0[1],I,J)*(f1[Ind2No(indPIp)]-f1[Ind2No(indPIm)])
	  *(f1[Ind2No(indPJp)]-f1[Ind2No(indPJm)])/(hPIm+hPIp)/(hPJm+hPJp);
      } // Cn
    }
  }

  uij /= coeff;

  return uij;
}

void JacobiPDE::PDE_solve(int maxstep, double tol, int n)
{
  double unext, u_norm;
  
  for (int step=0; step<maxstep; step++) {
    u_norm = 0;
    err = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:u_norm)
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
      cout << "\rerr for f1: " << setw(11) << left << err << "  step: " << step << flush;
    } else if (n == -2) {
      cout << "\rerr for g2: " << setw(11) << left << err << "  step: " << step << flush;
    }
    
    if (err < tol) {
      break;
    }
  }
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

double JacobiPDE::return_err()
{
  return err;
}

void JacobiPDE::export_fg(string filename)
{
  ofstream ofs(filename);

  vector<int> explist(volume,0), nextlist;
  double nextx, nextp, prex, prep;

  for (int number=0; number<volume; number++) {
    explist[number] = number;
  }

  
  for (int I=0; I<dim; I++) {
    nextx = site[0][I][0];
    nextp = site[1][I][0];
    prex = nextx;
    prep = nextp;
    
    for (int i=0; i<explist.size(); i++) {
      if (nextx > site[0][I][site[0][I].size()-1] && No2X(explist[i],I) < prex) {
	nextx = site[0][I][0];
      }
      
      if (No2X(explist[i],I) >= nextx || No2X(explist[i],I) == prex) {
	nextlist.push_back(explist[i]);
	prex = No2X(explist[i],I);

	if (nextx <= prex) {
	  nextx += exportdx[0][I];
	}
      }
    }
    explist.clear();
    explist = nextlist;
    nextlist.clear();

    for (int i=0; i<explist.size(); i++) {
      if (nextp > site[1][I][site[1][I].size()-1] && No2P(explist[i],I) < prep) {
	nextp = site[1][I][0];
      }
      
      if (No2P(explist[i],I) >= nextp || No2P(explist[i],I) == prep) {
	nextlist.push_back(explist[i]);
	prep = No2P(explist[i],I);

	if (nextp <= prep) {
	  nextp += exportdx[1][I];
	}
      }
    }
    explist.clear();
    explist = nextlist;
    nextlist.clear();
  }
  
  
  for (int i=0; i<explist.size(); i++) {
    for (int I=0; I<dim; I++) {
      ofs << No2X(explist[i],I) << ' ';
    }
    for (int I=0; I<dim; I++) {
      ofs << No2P(explist[i],I) << ' ';
    }
    ofs << f1[explist[i]] << ' '
	<< g2[explist[i]] << endl;
  }
}

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

#endif
