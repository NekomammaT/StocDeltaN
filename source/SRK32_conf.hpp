#ifndef INCLUDED_SRK32_hpp_
#define INCLUDED_SRK32_hpp_

#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include "MT.h"

using namespace std;

class SRKintegrater
{
protected:
  double t, t0,
    A0[3][3],A1[3][3],B0[3][3],B1[3][3],Alpha[3],Beta1[3],Beta2[3],C0[3],C1[3];
  vector<double> x,xi,dW,ax[3],H0x[3];
  vector< vector<double> > bx[3],Hkx[3],aIs,vIs;

public:
  SRKintegrater(){}
  SRKintegrater(vector<double> &Xi, double T0);
  void SRK2(double dt);
  void coeff(double dt, int step);
  double return_t();
  double return_phi(int I);
  double return_V();
  double vielbein(vector<double> &X, int I, int alpha);
  double eIsigma(vector<double> &X, int I);
  double eIs(vector<double> &X, int I, int alpha);
  virtual double V(vector<double> &X);
  virtual double VI(vector<double> &X, int I);
  virtual double metric(vector<double> &X, int I, int J);
  virtual double inversemetric(vector<double> &X, int I, int J);
  virtual double affine(vector<double> &X, int I, int J, int K);
  virtual double Dphi(vector<double> &X, int I);
  virtual double PhiNoise(vector<double> &X, int I, int alpha);
};

double Uniform();
double rand_normal(double mu, double sigma);



// --------------------- sample ---------------------------

double SRKintegrater::V(vector<double> &X) 
{
  double m1 = 0.01;
  double m2 = 0.1;

  return 1./2*m1*m1*X[0]*X[0] + 1./2*m2*m2*X[1]*X[1];
}

double SRKintegrater::VI(vector<double> &X, int I) 
{
  double m1 = 0.01;
  double m2 = 0.1;

  if (I == 0) {
    return m1*m1*X[0];
  } else {
    return m2*m2*X[1];
  }
}

double SRKintegrater::metric(vector<double> &X, int I, int J)
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

double SRKintegrater::inversemetric(vector<double> &X, int I, int J)
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

double SRKintegrater::affine(vector<double> &X, int I, int J, int K)
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

double SRKintegrater::Dphi(vector<double> &X, int I)
{
  double Dphi = 0;

  for (int J=0; J<X.size(); J++) {
    Dphi -= inversemetric(X,I,J)*VI(X,J)/V(X);
  }

  return Dphi;
}

double SRKintegrater::PhiNoise(vector<double> &X, int I, int alpha)
{
  return sqrt(V(X)/3.)/2./M_PI * vielbein(X,I,alpha);
}

// --------------------------------------------------------



SRKintegrater::SRKintegrater(vector<double> &Xi, double T0)
{
  init_genrand((unsigned)time(NULL));
  
  t = T0;
  x = Xi;
  t0 = T0;
  xi = Xi;

  for (int I=0; I<x.size(); I++) {
    aIs.push_back(x);
    vIs.push_back(x);
  }
  for (int I=0; I<x.size(); I++) {
    for (int alpha=0; alpha<x.size(); alpha++) {
      aIs[I][alpha] = rand_normal(0,1);
    }
  }
  
  for (int i=0; i<3; i++) {
    C0[i] = 0;
    C1[i] = 0;
    Alpha[i] = 0;
    Beta1[i] = 0;
    Beta2[i] = 0;
    ax[i] = x;
    H0x[i] = x;

    for (int j=0; j<3; j++) {
      A0[i][j] = 0;
      A1[i][j] = 0;
      B0[i][j] = 0;
      B1[i][j] = 0;
    }

    for (int alpha=0; alpha<x.size(); alpha++) {
      bx[i].push_back(x);
      Hkx[i].push_back(x);
    }
  }

  for (int alpha=0; alpha<x.size(); alpha++) {
    dW.push_back(0);
  }

  A0[1][0] = 1;
  A1[1][0] = 1;
  A1[2][0] = 1;
  B1[1][0] = 1;
  B1[2][0] = -1;
  C0[1] = 1;
  C1[1] = 1;
  C1[2] = 1;
  Alpha[0] = 1./2;
  Alpha[1] = 1./2;
  Beta1[0] = 1;
  Beta2[1] = 1./2;
  Beta2[2] = -1./2;
}

void SRKintegrater::SRK2(double dt)
{
  for (int alpha=0; alpha<dW.size(); alpha++) {
    dW[alpha] = rand_normal(0.,1.)*sqrt(dt);
  }

  for (int i=0; i<3; i++) {
    for (int I=0; I<x.size(); I++) {
      H0x[i][I] = x[I];

      for (int j=0; j<3; j++) {
	H0x[i][I] += A0[i][j]*ax[j][I]*dt;

	for (int alpha=0; alpha<dW.size(); alpha++) {
	  H0x[i][I] += B0[i][j]*bx[j][alpha][I]*dW[alpha];
	}
      }

      for (int beta=0; beta<dW.size(); beta++) {
	Hkx[i][beta][I] = x[I];

	for (int j=0; j<3; j++) {
	  Hkx[i][beta][I] += A1[i][j]*ax[j][I]*dt;

	  for (int alpha=0; alpha<dW.size(); alpha++) {
	    Hkx[i][beta][I] += B1[i][j]*bx[j][alpha][I]*dW[alpha]*dW[beta]/2./sqrt(dt);

	    if (alpha == beta) {
	      Hkx[i][beta][I] -= B1[i][j]*bx[j][alpha][I]*sqrt(dt)/2.;
	    }
	  }
	}
      }
    }

    coeff(dt,i);
  }

  for (int I=0; I<x.size(); I++) {
    for (int i=0; i<3; i++) {
      x[I] += Alpha[i]*ax[i][I]*dt;

      for (int alpha=0; alpha<dW.size(); alpha++) {
	x[I] += (Beta1[i]*dW[alpha]+Beta2[i]*sqrt(dt))*bx[i][alpha][I];
      }
    }
  }

  t += dt;
}

void SRKintegrater::coeff(double dt, int step)
{
  double T = t + C0[step]*dt;
  vector<double> X = H0x[step];

  for (int I=0; I<X.size(); I++) {
    ax[step][I] = Dphi(X,I);
  }


  T = t + C1[step]*dt;

  for (int alpha=0; alpha<dW.size(); alpha++) {
    X = Hkx[step][alpha];

    for (int I=0; I<X.size(); I++) {
      bx[step][alpha][I] = PhiNoise(X,I,alpha);
    }
  }
}

double SRKintegrater::return_t()
{
  return t;
}

double SRKintegrater::return_phi(int I)
{
  return x[I];
}

double SRKintegrater::return_V()
{
  return V(x);
}

double SRKintegrater::vielbein(vector<double> &X, int I, int alpha)
{
  if (alpha == 0) {
    return eIsigma(X,I);
  } else {
    return eIs(X,I,alpha);
  }
}

double SRKintegrater::eIsigma(vector<double> &X, int I)
{
  double NormN = 0;
  for (int J=0; J<X.size(); J++) {
    for (int K=0; K<X.size(); K++) {
      NormN += inversemetric(X,J,K)*VI(X,J)*VI(X,K);
    }
  }

  double eIsigma = 0;
  for (int J=0; J<X.size(); J++) {
    eIsigma += inversemetric(X,I,J)*VI(X,J);
  }

  return eIsigma/sqrt(NormN);
}

double SRKintegrater::eIs(vector<double> &X, int I, int alpha)
{
  vIs = aIs;
  
  double Norm;
  for (int al=0; al<alpha; al++) {
    Norm = 0;
    for (int J=0; J<X.size(); J++) {
      for (int K=0; K<X.size(); K++) {
	Norm += metric(X,J,K)*vielbein(X,J,al)*aIs[K][alpha];
      }
    }

    for (int J=0; J<X.size(); J++) {
      vIs[J][alpha] -= Norm*vielbein(X,J,al);
    }
  }

  Norm = 0;
  for (int J=0; J<X.size(); J++) {
    for (int K=0; K<X.size(); K++) {
      Norm += metric(X,J,K)*vIs[J][alpha]*vIs[K][alpha];
    }
  }

  return vIs[I][alpha]/sqrt(Norm);
}


double Uniform()
{
  return genrand_real3();
}

double rand_normal(double mu, double sigma)
{
  double z = sqrt(-2.*log(Uniform())) * sin(2.*M_PI*Uniform());
  return mu + sigma*z;
}


#endif
