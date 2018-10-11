#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "inflection_logp"

#define XMIN 0.3
#define XMAX 6
#define PMIN -20
#define PMAX -10
#define PC -12
#define HXOV 0.01
#define HXMIN (1e-4)
#define HXMAX 0.1
#define HPOV 0.01
#define HPMIN (1e-3)
#define HPMAX 0.1
#define MAXSTEP 100000
#define TOL 1e-10

#define LAMBDA 1
#define XI 2.3
#define PHIC 0.66
#define ALPHA (6*LAMBDA*PHIC/(3+XI*XI*PHIC*PHIC*PHIC*PHIC) - 4.3e-5)
#define M2 (LAMBDA*PHIC*PHIC*(3+XI*PHIC*PHIC)/(3+XI*XI*PHIC*PHIC*PHIC*PHIC))
#define V0 (1e-7)

#define HF (1.9e-5)
#define RHOC (3*HF*HF)

#define NOISEDIM 1
#define RECURSION 1000
#define PHIIN 5
#define PIN -12
#define TIMESTEP (1e-3)

#define DELTAN 0.01
#define NMAX 55


int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6; // start stop watch


  double h, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > sitepack[2];
  while (sitev <= XMAX) {
    h = min(max(fabs(sitev-PHIC)*HXOV, HXMIN), HXMAX);

    site.push_back(sitev);
    sitev += h;
  }
  sitepack[0].push_back(site);
  site.clear();

  sitev = PMIN;
  while (sitev <= PMAX) {
    h = min(max(fabs(sitev-PC)*HPOV, HPMIN), HPMAX);
    
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

  for (int I=0; I<dim; I++) {
    for (int J=0; J<dim; J++) {
      rho += 1./2*inversemetric(X,I,J)*exp(P[I]+P[J]);
    }
  }

  return sqrt(rho/3.);
}

double StocDeltaN::V(vector<double> &X)
{
  return V0/12.*(6*M2*X[0]*X[0] - 4*ALPHA*X[0]*X[0]*X[0] + 3*LAMBDA*X[0]*X[0]*X[0]*X[0])
    /(1+XI*X[0]*X[0])/(1+XI*X[0]*X[0]);
}

double StocDeltaN::VI(vector<double> &X, int I)
{
  return V0*X[0]*(M2*(3-3*XI*X[0]*X[0]) + X[0]*(3*LAMBDA*X[0]+ALPHA*(-3+XI*X[0]*X[0])))
    /3./(1+XI*X[0]*X[0])/(1+XI*X[0]*X[0])/(1+XI*X[0]*X[0]);
}

double StocDeltaN::VIJ(vector<double> &X, int I, int J)
{
  return V0*(3*M2*(1-8*XI*X[0]*X[0]+3*XI*XI*X[0]*X[0]*X[0]*X[0])
	     + X[0]*(-9*LAMBDA*X[0]*(-1+XI*X[0]*X[0])
		     - 2*ALPHA*(3-8*XI*X[0]*X[0]+XI*XI*X[0]*X[0]*X[0]*X[0])))
    /3./(1+XI*X[0]*X[0])/(1+XI*X[0]*X[0])/(1+XI*X[0]*X[0])/(1+XI*X[0]*X[0]);
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
      PiNoise += affine(X,K,I,J)*exp(P[K]-P[I])*vielbein(X,P,I,alpha) * H(X,P)/2./M_PI;
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

