#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "chaotic"

#define XMIN 0
#define XMAX 15
#define PMIN -16
#define PMAX -10
#define HX 0.1
#define HP 0.1
#define MAXSTEP 10000
#define TOL (1e-10)

#define MPHI (1e-5) //0.05

#define RHOC (MPHI*MPHI)

#define NOISEDIM 1
#define RECURSION 100
#define PHIIN 13
#define PIN -15
#define TIMESTEP (1e-2)

#define DELTAN 0.1
#define NMAX 40


int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6; // start stop watch


  double h = HX, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > sitepack[2];
  while (sitev <= XMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  sitepack[0].push_back(site);
  site.clear();

  sitev = PMIN;
  h = HP;
  while (sitev <= PMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  sitepack[1].push_back(site);
  site.clear();

  vector<double> xi = {PHIIN};
  vector<double> pi = {PIN};
  
  
  StocDeltaN sdn(MODEL,sitepack,RHOC,xi,pi,0,NOISEDIM,MAXSTEP,TOL,RECURSION,
		 TIMESTEP,NMAX,DELTAN);

  //sdn.sample();

  sdn.solve();

  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}






// --------------------- sample ---------------------------

double StocDeltaN::H(vector<double> &X, vector<double> &P)
{
  double rho = V(X);

  for (int I=0; I<X.size(); I++) {
    for (int J=0; J<X.size(); J++) {
      rho += 1./2*inversemetric(X,I,J)*exp(P[I]+P[J]);
    }
  }

  return sqrt(rho/3.);
}

double StocDeltaN::V(vector<double> &X)
{
  return 1./2*MPHI*MPHI*X[0]*X[0];
}

double StocDeltaN::VI(vector<double> &X, int I)
{
  return MPHI*MPHI*X[0];
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

double StocDeltaN::PhiNoise(vector<double> &X, vector<double> &P, int I, int alpha)
{
  return H(X,P)/2./M_PI * vielbein(X,P,I,alpha);
}

double StocDeltaN::PiNoise(vector<double> &X, vector<double> &P, int I, int alpha)
{
  double PiNoise = 0;

  for (int J=0; J<dim; J++) {
    for (int K=0; K<dim; K++) {
      PiNoise += affine(X,K,I,J)*exp(P[K]-P[I])*PhiNoise(X,P,J,alpha);
    }
  }

  return PiNoise;
}

double StocDeltaN::Dphi(vector<double> &X, vector<double> &P, int I)
{
  double Dphi = 0;

  for (int J=0; J<dim; J++) {
    Dphi -= inversemetric(X,I,J)*exp(P[J])/H(X,P);
  }

  return Dphi;
}

double StocDeltaN::Dpi(vector<double> &X, vector<double> &P, int I)
{
  double Hubble = H(X,P);
  double Dpi = -3+VI(X,I)*exp(-P[I])/Hubble;

  for (int J=0; J<dim; J++) {
    for (int K=0; K<dim; K++) {
      for (int L=0; L<dim; L++) {
	Dpi -= affine(X,K,I,J)*inversemetric(X,J,L)*exp(P[K]+P[L]-P[I])/Hubble;
      }
    }
  }

  return Dpi;
}

double StocDeltaN::Dphiphi(vector<double> &X, vector<double> &P, int I, int J)
{
  return H(X,P)*H(X,P)/4./M_PI/M_PI * inversemetric(X,I,J);
}

double StocDeltaN::Dphipi(vector<double> &X, vector<double> &P, int I, int J)
{
  double Dphipi = 0;

  for (int K=0; K<dim; K++) {
    for (int L=0; L<dim; L++) {
      Dphipi += affine(X,K,J,L)*exp(P[K]-P[J])*Dphiphi(X,P,I,L);
    }
  }

  return Dphipi;
}

double StocDeltaN::Dpipi(vector<double> &X, vector<double> &P, int I, int J)
{
  double Dpipi = 0;

  for (int K=0; K<dim; K++) {
    for (int L=0; L<dim; L++) {
      for (int M=0; M<dim; M++) {
	for (int N=0; N<dim; N++) {
	  Dpipi += affine(X,K,I,L)*affine(X,M,J,N)*exp(P[K]+P[M]-P[I]-P[J])*Dphiphi(X,P,L,N);
	}
      }
    }
  }

  return Dpipi;
}

// --------------------------------------------------------
