#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "chaotic_logp" // model name, momentum = log(-pi)

// ---------- box size & step h ------------
#define XMIN 0
#define XMAX 15
#define PMIN -16
#define PMAX -10
#define HX 0.1
#define HP 0.1
// -----------------------------------------

// ---------- for PDE ----------------------
#define MAXSTEP 10000 // max recursion
#define TOL (1e-10) // tolerance
// -----------------------------------------

#define MPHI (1e-5) // potential parameter

#define RHOC (MPHI*MPHI) // end of inflation

// ---------- for SDE ----------------------
#define NOISEDIM 1 // d.o.f. of noise
#define RECURSION 100 // recursion for power spectrum
#define PHIIN 13 // i.c. for phi
#define PIN -15 // i.c. for log(-pi)
#define TIMESTEP (1e-2) // time step: delta N
// -----------------------------------------

// ---------- for power spectrum -----------
#define DELTAN 0.1 // calc. PS every DELTAN e-folds
#define NMAX 40 // calc. PS for 0--NMAX e-folds
// -----------------------------------------


int main(int argc, char** argv)
{
  // ---------- start stopwatch --------------
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  // -----------------------------------------

  // ---------- set box step h ---------------
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
  // ------------------------------------------

  // ---------- set i.c. for sample paths -----
  vector<double> xi = {PHIIN};
  vector<double> pi = {PIN};
  // ------------------------------------------
  
  StocDeltaN sdn(MODEL,sitepack,RHOC,xi,pi,0,NOISEDIM,MAXSTEP,TOL,RECURSION,
		 TIMESTEP,NMAX,DELTAN); // declare the system

  //sdn.sample(); // show 1 sample path

  sdn.solve(); // solve PDE & SDE and obtain power spectrum

  // --------- stop stopwatch ----------------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // -----------------------------------------
}


// ---------- Lagrangian params. and diff. coeff.  X[0]=phi, P[0]=log(-pi) ----------

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

double StocDeltaN::VI(vector<double> &X, int I) // \partial_I V
{
  return MPHI*MPHI*X[0];
}

double StocDeltaN::metric(vector<double> &X, int I, int J) // G_IJ
{
  return 1;
}

double StocDeltaN::inversemetric(vector<double> &X, int I, int J) // G^IJ
{
  return 1;
}

double StocDeltaN::affine(vector<double> &X, int I, int J, int K) // \Gamma^I_JK
{
  return 0;
}

double StocDeltaN::PhiNoise(vector<double> &X, vector<double> &P, int I, int alpha) // gQ^I_a
{
  return H(X,P)/2./M_PI * vielbein(X,P,I,alpha);
}

double StocDeltaN::PiNoise(vector<double> &X, vector<double> &P, int I, int alpha) // gP_Ia
{
  double PiNoise = 0;

  for (int J=0; J<dim; J++) {
    for (int K=0; K<dim; K++) {
      PiNoise += affine(X,K,I,J)*exp(P[K]-P[I])*PhiNoise(X,P,J,alpha);
    }
  }

  return PiNoise;
}

double StocDeltaN::Dphi(vector<double> &X, vector<double> &P, int I) // Dphi^I
{
  double Dphi = 0;

  for (int J=0; J<dim; J++) {
    Dphi -= inversemetric(X,I,J)*exp(P[J])/H(X,P);
  }

  return Dphi;
}

double StocDeltaN::Dpi(vector<double> &X, vector<double> &P, int I) // Dpi_I
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

double StocDeltaN::Dphiphi(vector<double> &X, vector<double> &P, int I, int J) // Dphiphi^IJ
{
  return H(X,P)*H(X,P)/4./M_PI/M_PI * inversemetric(X,I,J);
}

double StocDeltaN::Dphipi(vector<double> &X, vector<double> &P, int I, int J) // Dphipi^I_J
{
  double Dphipi = 0;

  for (int K=0; K<dim; K++) {
    for (int L=0; L<dim; L++) {
      Dphipi += affine(X,K,J,L)*exp(P[K]-P[J])*Dphiphi(X,P,I,L);
    }
  }

  return Dphipi;
}

double StocDeltaN::Dpipi(vector<double> &X, vector<double> &P, int I, int J) // Dpipi_IJ
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
