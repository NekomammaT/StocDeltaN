#include "SRK32.hpp"
#include "MT.h"

// --------------------- sample ---------------------------

double SRKintegrater::H(vector<double> &X, vector<double> &P)
{
  double rho = V(X);

  for (int I=0; I<X.size(); I++) {
    for (int J=0; J<X.size(); J++) {
      rho += 1./2*inversemetric(X,I,J)*P[I]*P[J];
    }
  }

  return sqrt(rho/3.);
}

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

double SRKintegrater::Dphi(vector<double> &X, vector<double> &P, int I)
{
  double Dphi = 0;

  for (int J=0; J<X.size(); J++) {
    Dphi += inversemetric(X,I,J)*P[J]/H(X,P);
  }

  return Dphi;
}

double SRKintegrater::Dpi(vector<double> &X, vector<double> &P, int I)
{
  double Hubble = H(X,P);
  double Dpi = -3*P[I]-VI(X,I)/Hubble;

  for (int J=0; J<X.size(); J++) {
    for (int K=0; K<X.size(); K++) {
      for (int L=0; L<X.size(); L++) {
	Dpi += affine(X,K,I,J)*inversemetric(X,J,L)*P[K]*P[L]/Hubble;
      }
    }
  }

  return Dpi;
}

double SRKintegrater::PhiNoise(vector<double> &X, vector<double> &P, int I, int alpha)
{
  return H(X,P)/2./M_PI * vielbein(X,P,I,alpha);
}

double SRKintegrater::PiNoise(vector<double> &X, vector<double> &P, int I, int alpha)
{
  double PiNoise = 0;

  for (int J=0; J<X.size(); J++) {
    for (int K=0; K<X.size(); K++) {
      PiNoise += affine(X,K,I,J)*P[K]*vielbein(X,P,I,alpha) * H(X,P)/2./M_PI;
    }
  }

  return PiNoise;
}


// --------------------------------------------------------



SRKintegrater::SRKintegrater(vector<double> &Xi, vector<double> &Pi, double T0,
			     int NoiseDim)
{
  init_genrand((unsigned)time(NULL));
  
  t = T0;
  x = Xi;
  p = Pi;
  t0 = T0;
  xi = Xi;
  pi = Pi;

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
    ap[i] = p;
    H0x[i] = x;
    H0p[i] = p;

    for (int j=0; j<3; j++) {
      A0[i][j] = 0;
      A1[i][j] = 0;
      B0[i][j] = 0;
      B1[i][j] = 0;
    }

    for (int alpha=0; alpha<NoiseDim; alpha++) {
      bx[i].push_back(x);
      bp[i].push_back(p);
      Hkx[i].push_back(x);
      Hkp[i].push_back(p);
    }
  }

  for (int alpha=0; alpha<NoiseDim; alpha++) {
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
      H0p[i][I] = p[I];

      for (int j=0; j<3; j++) {
	H0x[i][I] += A0[i][j]*ax[j][I]*dt;
	H0p[i][I] += A0[i][j]*ap[j][I]*dt;

	for (int alpha=0; alpha<dW.size(); alpha++) {
	  H0x[i][I] += B0[i][j]*bx[j][alpha][I]*dW[alpha];
	  H0p[i][I] += B0[i][j]*bp[j][alpha][I]*dW[alpha];
	}
      }

      for (int beta=0; beta<dW.size(); beta++) {
	Hkx[i][beta][I] = x[I];
	Hkp[i][beta][I] = p[I];

	for (int j=0; j<3; j++) {
	  Hkx[i][beta][I] += A1[i][j]*ax[j][I]*dt;
	  Hkp[i][beta][I] += A1[i][j]*ap[j][I]*dt;

	  for (int alpha=0; alpha<dW.size(); alpha++) {
	    Hkx[i][beta][I] += B1[i][j]*bx[j][alpha][I]*dW[alpha]*dW[beta]/2./sqrt(dt);
	    Hkp[i][beta][I] += B1[i][j]*bp[j][alpha][I]*dW[alpha]*dW[beta]/2./sqrt(dt);

	    if (alpha == beta) {
	      Hkx[i][beta][I] -= B1[i][j]*bx[j][alpha][I]*sqrt(dt)/2.;
	      Hkp[i][beta][I] -= B1[i][j]*bp[j][alpha][I]*sqrt(dt)/2.;
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
      p[I] += Alpha[i]*ap[i][I]*dt;

      for (int alpha=0; alpha<dW.size(); alpha++) {
	x[I] += (Beta1[i]*dW[alpha]+Beta2[i]*sqrt(dt))*bx[i][alpha][I];
	p[I] += (Beta1[i]*dW[alpha]+Beta2[i]*sqrt(dt))*bp[i][alpha][I];
      }
    }
  }

  t += dt;
}

void SRKintegrater::coeff(double dt, int step)
{
  double T = t + C0[step]*dt;
  vector<double> X = H0x[step];
  vector<double> P = H0p[step];
  double Hubble = H(X,P);

  for (int I=0; I<X.size(); I++) {
    ax[step][I] = Dphi(X,P,I);
    ap[step][I] = Dpi(X,P,I);
  }


  T = t + C1[step]*dt;

  for (int alpha=0; alpha<dW.size(); alpha++) {
    X = Hkx[step][alpha];
    P = Hkp[step][alpha];

    for (int I=0; I<X.size(); I++) {
      bx[step][alpha][I] = PhiNoise(X,P,I,alpha);
      bp[step][alpha][I] = PiNoise(X,P,I,alpha);
    }
  }
}


double SRKintegrater::e1(vector<double> &X, vector<double> &P)
{
  double e1 = 0;
  for (int I=0; I<P.size(); I++) {
    for (int J=0; J<P.size(); J++) {
      e1 += inversemetric(X,I,J)*P[I]*P[J];
    }
  }
  e1 /= 2*H(X,P)*H(X,P);

  return e1;
}

double SRKintegrater::return_H()
{
  return H(x,p);
}

double SRKintegrater::return_t()
{
  return t;
}

double SRKintegrater::return_phi(int I)
{
  return x[I];
}

double SRKintegrater::return_pi(int I)
{
  return p[I];
}

double SRKintegrater::return_e1()
{
  return e1(x,p);
}

double SRKintegrater::vielbein(vector<double> &X, vector<double> &P, int I, int alpha)
{
  if (alpha == 0) {
    return eIsigma(X,P,I);
  } else {
    return eIs(X,P,I,alpha);
  }
}

double SRKintegrater::eIsigma(vector<double> &X, vector<double> &P, int I)
{
  double NormN = 0;
  for (int J=0; J<X.size(); J++) {
    for (int K=0; K<X.size(); K++) {
      NormN += inversemetric(X,J,K)*P[J]*P[K];
    }
  }

  double eIsigma = 0;
  for (int J=0; J<X.size(); J++) {
    eIsigma += inversemetric(X,I,J)*P[J];
  }

  return eIsigma/sqrt(NormN);
}

double SRKintegrater::eIs(vector<double> &X, vector<double> &P, int I, int alpha)
{
  vIs = aIs;
  
  double Norm;
  for (int al=0; al<alpha; al++) {
    Norm = 0;
    for (int J=0; J<X.size(); J++) {
      for (int K=0; K<X.size(); K++) {
	Norm += metric(X,J,K)*vielbein(X,P,J,al)*aIs[K][alpha];
      }
    }

    for (int J=0; J<X.size(); J++) {
      vIs[J][alpha] -= Norm*vielbein(X,P,J,al);
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
