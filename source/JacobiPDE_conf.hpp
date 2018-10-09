#ifndef INCLUDED_JacobiPDE_hpp_
#define INCLUDED_JacobiPDE_hpp_

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

class JacobiPDE
{
protected:
  vector<double> f1, f1_next, g2, g2_next; // argument: number
  vector< vector<double> > site, hI; // argument: I, siteNo
  vector<double> FPoint; // argument: I
  vector<bool> Omega; // argument: number
  vector<int> ind; // argument: I
  vector<int> siteNo; // argument: I
  double rhoc, err;
  int dim, volume, pvol;
  
public:
  JacobiPDE(){}
  JacobiPDE(vector< vector<double> > &Site, double Rhoc);
  double PDE(int number, int n);
  void PDE_solve(int maxstep, double tol, int n);
  int Ind2No(vector<int> &index);
  int No2XInd(int number, int I);
  double No2X(int number, int I);
  int ceilX(vector<double> &X, int I);
  double Interpolation_f(vector<double> &X, vector<double> &f);
  double return_f1(vector<int> &index);
  double return_g2(vector<int> &index);
  double return_err();
  void export_fg(string filename);
  virtual double V(vector<double> &X);
  virtual double VI(vector<double> &X, int I);
  virtual double metric(vector<double> &X, int I, int J);
  virtual double inversemetric(vector<double> &X, int I, int J);
  virtual double affine(vector<double> &X, int I, int J, int K);
  virtual double Dphi(vector<double> &X, int I);
  virtual double Dphiphi(vector<double> &X, int I, int J);
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

double JacobiPDE::Dphi(vector<double> &X, int I)
{
  double Dphi = 0;

  for (int J=0; J<dim; J++) {
    Dphi -= inversemetric(X,I,J)*VI(X,J)/V(X);
  }

  return Dphi;
}

double JacobiPDE::Dphiphi(vector<double> &X, int I, int J)
{
  return V(X)/12./M_PI/M_PI * inversemetric(X,I,J);
}

// -------------------------------------------------------



JacobiPDE::JacobiPDE(vector< vector<double> > &Site, double Rhoc)
{
  srand((unsigned)time(NULL));

#ifdef _OPENMP
  cout<<"OpenMP : Enabled (Max # of threads = "<<omp_get_max_threads()<<")"<<endl;
#endif

  rhoc = Rhoc;
  site = Site;
  hI = site;
  dim = site.size();
  siteNo = vector<int>(dim,0);
  ind = vector<int>(dim,0);
  FPoint = vector<double>(dim,0);

  for (int I=0; I<dim; I++) {
    siteNo[I] = site[I].size();
    
    for (int s=0; s<site[I].size()-1; s++) {
      hI[I][s] = site[I][s+1] - site[I][s];
    }
  }

  volume = 1;
  for (int I=0; I<dim; I++) {
    volume *= siteNo[I];
  }

  f1 = vector<double>(volume,0);
  f1_next = f1;
  g2 = f1;
  g2_next = f1;
  Omega = vector<bool>(volume,true);

  for (int number=0; number<volume; number++) {
    for (int I=0; I<dim; I++) {
      FPoint[I] = No2X(number,I);
    }
    
    if (V(FPoint) < rhoc) {
      Omega[number] = false;
      f1[number] = 0;
      g2[number] = 0;
    } else {
      Omega[number] = true;
      f1[number] = rand()%10;
      g2[number] = (rand()%10)/10.;
    }
  }
}

double JacobiPDE::PDE(int number, int n)
{
  vector<int> ind0;
  vector<double> FPoint0;

  ind0 = ind;
  FPoint0 = FPoint;
  
  for (int I=0; I<dim; I++) {
    ind0[I] = No2XInd(number,I);
    FPoint0[I] = No2X(number,I);
  }
  
  vector<int> indXIm, indXIp, indXJm, indXJp, indXXmm, indXXmp, indXXpm, indXXpp;
  double uXIm, uXIp, uXXmm, uXXmp, uXXpm, uXXpp, hXIp, hXIm, hXJp, hXJm;
  double uij = 0, coeff = 0;

  if (n == 1) {
    uij -= 1;
  } // Cn
  
  for (int I=0; I<dim; I++) {
    indXIm = ind0;
    indXIp = ind0;
    
    // -------------- reflecting b.c. -----------------
    
    if (indXIm[I] == 0) {
      indXIm[I]++;
      hXIm = hI[I][ind0[I]];
    } else {
      indXIm[I]--;
      hXIm = hI[I][ind0[I]-1];
    }

    if (indXIp[I] == siteNo[I]-1) {
      indXIp[I]--;
      hXIp = hI[I][ind0[I]-1];
    } else {
      indXIp[I]++;
      hXIp = hI[I][ind0[I]];
    }

    // ----------------------------------------------

    if (n == 1) {
      uXIm = f1[Ind2No(indXIm)];
      uXIp = f1[Ind2No(indXIp)];
    } else if (n == -2) {
      uXIm = g2[Ind2No(indXIm)];
      uXIp = g2[Ind2No(indXIp)];
    }

    if (Dphi(FPoint0,I) < 0 ) {
      uij += Dphi(FPoint0,I)/hXIm * uXIm;
      coeff += Dphi(FPoint0,I)/hXIm;
    } else {
      uij -= Dphi(FPoint0,I)/hXIp * uXIp;
      coeff -= Dphi(FPoint0,I)/hXIp;
    }

    uij -= Dphiphi(FPoint0,I,I)*(uXIp*hXIm+uXIm*hXIp)/hXIp/hXIm/(hXIm+hXIp);
    coeff -= Dphiphi(FPoint0,I,I)/hXIp/hXIm;
    
    
    for (int J=0; J<dim; J++) {
      indXJm = ind0;
      indXJp = ind0;
      indXXmm = ind0;
      indXXmp = ind0;
      indXXpm = ind0;
      indXXpp = ind0;
      
      // -------------------- reflecting b.c. --------------------
      
      if (indXJm[J] == 0) {
	indXJm[J]++;
	hXJm = hI[J][ind0[J]];
      } else {
	indXJm[J]--;
	hXJm = hI[J][ind0[J]-1];
      }
      if (indXJp[J] == siteNo[J]-1) {
	indXJp[J]--;
	hXJp = hI[J][ind0[J]-1];
      } else {
	indXJp[J]++;
	hXJp = hI[J][ind0[J]];
      }
      
      if (indXXmm[I] == 0) {
	indXXmm[I]++;
      } else {
	indXXmm[I]--;
      }
      if (indXXmp[I] == 0) {
	indXXmp[I]++;
      } else {
	indXXmp[I]--;
      }
      if (indXXpm[I] == siteNo[I]-1) {
	indXXpm[I]--;
      } else {
	indXXpm[I]++;
      }
      if (indXXpp[I] == siteNo[I]-1) {
	indXXpp[I]--;
      } else {
	indXXpp[I]++;
      }
      
      if (indXXmm[J] == 0) {
	indXXmm[J]++;
      } else {
	indXXmm[J]--;
      }
      if (indXXpm[J] == 0) {
	indXXpm[J]++;
      } else {
	indXXpm[J]--;
      }
      if (indXXmp[J] == siteNo[J]-1) {
	indXXmp[J]--;
      } else {
	indXXmp[J]++;
      }
      if (indXXpp[J] == siteNo[J]-1) {
	indXXpp[J]--;
      } else {
	indXXpp[J]++;
      }

      // ----------------------------------

      if (n == 1) {
	uXXmm = f1[Ind2No(indXXmm)];
	uXXmp = f1[Ind2No(indXXmp)];
	uXXpm = f1[Ind2No(indXXpm)];
	uXXpp = f1[Ind2No(indXXpp)];
      } else if (n == -2) {
	uXXmm = g2[Ind2No(indXXmm)];
	uXXmp = g2[Ind2No(indXXmp)];
	uXXpm = g2[Ind2No(indXXpm)];
	uXXpp = g2[Ind2No(indXXpp)];
      }
             
      if (J != I) {
	uij -= 1./2*Dphiphi(FPoint0,I,J)*(uXXpp-uXXpm-uXXmp+uXXmm)/(hXIm+hXIp)/(hXJm+hXJp);
      }

      if (n == -2) {
	uij -= Dphiphi(FPoint0,I,J)*(f1[Ind2No(indXIp)]-f1[Ind2No(indXIm)])
	  *(f1[Ind2No(indXJp)]-f1[Ind2No(indXJm)])/(hXIm+hXIp)/(hXJm+hXJp);
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

int JacobiPDE::Ind2No(vector<int> &index)
{
  int number = 0;
  int xtemp;

  for (int I=0; I<dim; I++) {
    xtemp = index[dim-1-I];
    
    for (int J=0; J<I; J++) {
      xtemp *= siteNo[dim-1-J];
    }
    
    number += xtemp;
  }
  
  return number;
}

int JacobiPDE::No2XInd(int number, int I)
{
  int index = number;

  for (int J=dim-1; J>I; J--) {
    index /= siteNo[J];
  }

  if (I > 0) {
    index %= siteNo[I];
  }

  return index;
}

double JacobiPDE::No2X(int number, int I)
{
  return site[I][No2XInd(number,I)];
}

int JacobiPDE::ceilX(vector<double> &X, int I)
{
  int xind = 0;

  while (site[I][xind] <= X[I]) {
    xind++;
  }

  return xind;
}

double JacobiPDE::Interpolation_f(vector<double> &X, vector<double> &f) 
{
  int pmcheckX;
  double weight, intf = 0;
  vector< vector<int> > near_site;
  vector<double> fdata;
  
  for (int xnum=0; xnum<pow(2,dim); xnum++) {
    weight = 1;
    
    for (int I=0; I<dim; I++) {
      pmcheckX = xnum;
      
      for (int J=0; J<I; J++) {
	pmcheckX /= 2;
      }
	
      if (I < dim-1) {
	pmcheckX %= 2;
      }

      if (pmcheckX == 1) {
	ind[I] = ceilX(X,I);
	weight *= (X[I]-site[I][ceilX(X,I)-1])/hI[I][ceilX(X,I)-1];
      } else {
	ind[I] = ceilX(X,I)-1;
	weight *= (site[I][ceilX(X,I)]-X[I])/hI[I][ceilX(X,I)-1];
      }
    }
    
    intf += f[Ind2No(ind)]*weight;
  }

  return intf;
}

double JacobiPDE::return_f1(vector<int> &index)
{
  return f1[Ind2No(index)];
}

double JacobiPDE::return_g2(vector<int> &index)
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

  for (int number=0; number<volume; number++) {
    for (int I=0; I<dim; I++) {
      ofs << No2X(number,I) << ' ';
    }
    ofs << f1[number] << ' '
	<< g2[number] << endl;
  }
}

#endif
