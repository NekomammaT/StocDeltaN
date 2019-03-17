#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "USR"

#define XMIN 4.23e-2
#define XMAX 0.2
#define PMIN -(2e-6)
#define PMAX 0
#define HXMIN (1e-7)
#define HXOV (1./20)
#define HPMIN (1e-11)
#define HPOV (1./20)

#define MAXSTEP 100000
#define TOL 1e-10

#define V0 (1e-10)
#define XF (4.232e-2)

#define RHOC V0

#define NOISEDIM 1
#define RECURSION 100
#define XIN 0.1 
#define PIN -(1e-6)
#define TIMESTEP (1e-4)

#define DELTAN 0.01
#define NMAX 3


int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;

  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;

  double h, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > sitepack[2];
  while (sitev <= XMAX) {
    h = max(fabs(sitev-XF)*HXOV,HXMIN);

    site.push_back(sitev);
    sitev += h;
  }
  sitepack[0].push_back(site);
  site.clear();

  sitev = PMIN;
  while (sitev <= PMAX) {
    h = max(fabs(sitev)*HPOV,HPMIN);

    site.push_back(sitev);
    sitev += h;
  }
  sitepack[1].push_back(site);
  site.clear();

  vector<double> xi = {XIN};
  vector<double> pi = {PIN};

  StocDeltaN sdn(MODEL,sitepack,RHOC,xi,pi,0,NOISEDIM,MAXSTEP,TOL,RECURSION,
		 TIMESTEP,NMAX,DELTAN);

  //sdn.sample();
  //sdn.sample_loglogplot();

  sdn.solve();
  sdn.calP_plot();

  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
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
  if (X[0] > XF) {
    return V0;
  } else {
    return 0;
  }
}

double StocDeltaN::VI(vector<double> &X, int I)
{
  return 0;
}

double StocDeltaN::metric(vector<double> &X, int I, int J)
{
  return 1;
}

double StocDeltaN::inversemetric(vector<double> &X, int I, int J)
{
  return 1;
}

double StocDeltaN::affine(vector<double> &X, int I, int J, int K)
{
  return 0;
}

double StocDeltaN::derGamma(vector<double> &X, int I, int J, int K, int L)
{
  return 0;
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

  for (int J=0; J<X.size(); J++) {
    for (int K=0; K<X.size(); K++) {
      Dpi += affine(X,J,I,K)*Dphipi(X,P,K,J);
      
      for (int S=0; S<X.size(); S++) {
	Dpi += affine(X,S,I,J)*inversemetric(X,J,K)*P[K]*P[S]/Hubble
	  +1./2*derGamma(X,S,I,J,K)*P[S]*Dphiphi(X,P,J,K);

	for (int R=0; R<X.size(); R++) {
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
