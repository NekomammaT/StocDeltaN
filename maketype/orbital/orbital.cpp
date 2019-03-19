#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "orbital"

#define RHOMIN 0.99
#define RHOMAX 1.01
#define HRHO (1e-3)
#define THETAMIN 1
#define THETAMAX 13
#define HTHETA 1
#define RPIMIN -(4e-11)
#define RPIMAX (1e-10)
#define HRPI (1e-11)
#define TPIMIN -(9e-6)
#define TPIMAX -(7e-6)
#define HTPI (1e-6)

#define MAXSTEP 100000
#define TOL 1e-10

#define MM (1e-5)

#define RHOC (MM*MM)

#define NOISEDIM 2
#define RECURSION 100
#define RHOIN 1
#define RPIIN 0
#define THETAIN 10
#define TPIIN -sqrt(2./3)*MM
#define TIMESTEP (1e-2)

#define DELTAN 0.1
#define NMAX 24


int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6; // start stop watch


  double h = HRHO, sitev = RHOMIN;
  vector<double> site;
  vector< vector<double> > sitepack[2];
  while (sitev <= RHOMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  sitepack[0].push_back(site);
  site.clear();

  h = HTHETA, sitev = THETAMIN;
  while (sitev <= THETAMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  sitepack[0].push_back(site);
  site.clear();

  h = HRPI, sitev = RPIMIN;
  while (sitev <= RPIMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  sitepack[1].push_back(site);
  site.clear();

  h = HTPI, sitev = TPIMIN;
  while (sitev <= TPIMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  sitepack[1].push_back(site);
  site.clear();

  vector<double> xi = {RHOIN,THETAIN};
  vector<double> pi = {RPIIN,TPIIN};

  StocDeltaN sdn(MODEL,sitepack,RHOC,xi,pi,0,NOISEDIM,MAXSTEP,TOL,RECURSION,
		 TIMESTEP,NMAX,DELTAN);

  //sdn.sample();
  //sdn.sample_plot();

  /*
  vector<double> XP[2];
  XP[0] = vector<double> {0.991,7};
  XP[1] = vector<double> {5e-11,-7e-6};
  cout << sdn.Dphi(XP[0],XP[1],0) << endl;
  */
  
  sdn.solve();
  //sdn.calP_plot();
  
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}


double StocDeltaN::H(vector<double> &X, vector<double> &P)
{
  double rho = V(X);

  for (int I=0; I<X.size(); I++) {
    for (int J=0; J<X.size(); J++) {
      rho += 1./2*inversemetric(X,I,J)*P[I]*P[J];
    }
  }

  return sqrt(rho/3.);
}

double StocDeltaN::V(vector<double> &X)
{
  return 1./2*MM*MM*(X[1]*X[1]-2./3/X[0]/X[0]);
}

double StocDeltaN::VI(vector<double> &X, int I) // \partial_I V
{
  if (I == 0) {
    return 2*MM*MM/3./X[0]/X[0];
  } else {
    return MM*MM*X[1];
  }
}

double StocDeltaN::metric(vector<double> &X, int I, int J) // G_IJ
{
  if (I == 0 && J == 0) {
    return 1;
  } else if (I == 1 && J == 1) {
    return X[0]*X[0];
  } else {
    return 0;
  }
}

double StocDeltaN::inversemetric(vector<double> &X, int I, int J) // G^IJ
{
  if (I == 0 && J == 0) {
    return 1;
  } else if (I == 1 && J == 1) {
    return 1./X[0]/X[0];
  } else {
    return 0;
  }
}

double StocDeltaN::affine(vector<double> &X, int I, int J, int K) // \Gamma^I_JK
{
  if (I == 0 && J == 1 && K == 1) {
    return -X[0];
  } else if (I == 1 && J != K) {
    return 1./X[0];
  } else {
    return 0;
  }
}

double StocDeltaN::derGamma(vector<double> &X, int I, int J, int K, int L)
{
  if (L == 0) {
    if (I == 0 && J == 1 && K == 1) {
      return -1;
    } else if (I == 1 && J != K) {
      return -1./X[0]/X[0];
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

double StocDeltaN::PhiNoise(vector<double> &X, vector<double> &P, int I, int alpha)
{
  return H(X,P)/2./M_PI * vielbein(X,P,I,alpha);
}

double StocDeltaN::PiNoise(vector<double> &X, vector<double> &P, int I, int alpha)
{
  double PiNoise = 0;

  for (int J=0; J<dim; J++) {
    for (int K=0; K<dim; K++) {
      PiNoise += affine(X,K,I,J)*P[K]*PhiNoise(X,P,J,alpha);
    }
  }

  return PiNoise;
}

double StocDeltaN::Dphi(vector<double> &X, vector<double> &P, int I)
{
  double Dphi = 0;

  for (int J=0; J<dim; J++) {
    Dphi += inversemetric(X,I,J)*P[J]/H(X,P);

    for (int K=0; K<dim; K++) {
      Dphi -= 1./2*affine(X,I,J,K)*Dphiphi(X,P,J,K);
    }
  }

  return Dphi;
}

double StocDeltaN::Dpi(vector<double> &X, vector<double> &P, int I)
{
  double Hubble = H(X,P);
  double Dpi = -3*P[I]-VI(X,I)/Hubble;

  for (int J=0; J<dim; J++) {
    for (int K=0; K<dim; K++) {
      Dpi += affine(X,J,I,K)*Dphipi(X,P,K,J);

      for (int S=0; S<dim; S++) {
	Dpi += affine(X,S,I,J)*inversemetric(X,J,K)*P[K]*P[S]/Hubble
	  +1./2*derGamma(X,S,I,J,K)*P[S]*Dphiphi(X,P,J,K);

	for (int R=0; R<dim; R++) {
	  Dpi -= 1./2*(affine(X,R,J,K)*affine(X,S,I,R) + affine(X,R,I,J)*affine(X,S,K,R))
	    *P[S]*Dphiphi(X,P,J,K);
	}
      }
    }
  }

  return Dpi;
}

double StocDeltaN::Dphiphi(vector<double> &X, vector<double> &P, int I, int J)
{
  double Dphiphi = 0;

  for (int alpha=0; alpha<dW.size(); alpha++) {
    Dphiphi += PhiNoise(X,P,I,alpha)*PhiNoise(X,P,J,alpha);
  }

  return Dphiphi;
}

double StocDeltaN::Dphipi(vector<double> &X, vector<double> &P, int I, int J)
{
  double Dphipi = 0;

  for (int alpha=0; alpha<dW.size(); alpha++) {
    Dphipi += PhiNoise(X,P,I,alpha)*PiNoise(X,P,J,alpha);
  }

  return Dphipi;
}

double StocDeltaN::Dpipi(vector<double> &X, vector<double> &P, int I, int J)
{
  double Dpipi = 0;

  for (int alpha=0; alpha<dW.size(); alpha++) {
    Dpipi += PiNoise(X,P,I,alpha)*PiNoise(X,P,J,alpha);
  }

  return Dpipi;
}



