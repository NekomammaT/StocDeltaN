#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "hybrid_SR" // model name

// ---------- box size & step h ------------
#define PHIMIN 0.1409
#define PHIMAX 0.142
#define PSIMIN -(1e-3)
#define PSIMAX (1e-3)
#define HPHI (1e-5)
#define HPSIOPSI (1e-2) // hpsi/|psi|
#define HPSIMIN (1e-10)
// -----------------------------------------

// ---------- for PDE ----------------------
#define MAXSTEP 100000 // max recursion
#define TOL 1e-10 // tolerance
// -----------------------------------------

// ---------- potential parameter ----------
#define AS (2.189e-9)
#define PI2 50
#define MM 0.1
#define PHIC (0.1*sqrt(2))
#define MU2 11
#define MU1 (PI2/MM/MM/PHIC)
#define LAMBDA4 (AS*12*M_PI*M_PI/MU1/MU1)
// -----------------------------------------

#define RHOC (2.074038e-16) // end of inflation

// ---------- for SDE ----------------------
#define RECURSION 10000 // recursion for power spectrum
#define PHIIN 0.1418 // i.c. for phi
#define PSIIN 0 // i.c. for psi
#define TIMESTEP (1e-3) // time step : delta N
// -----------------------------------------

// ---------- for power spectrum -----------
#define DELTAN 0.1 // calc. PS every DELTAN e-folds
#define NMAX 28 // calc. PS for 0--NMAX e-folds
// -----------------------------------------


int main(int argc, char** argv)
{
  // ---------- start stop watch ----------
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  // --------------------------------------


  // ---------- set box ---------------
  double h = HPHI, sitev = PHIMIN;
  vector<double> site;
  vector< vector<double> > xsite;
  vector< vector< vector<double> > > sitepack;
  while (sitev <= PHIMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();

  sitev = PSIMIN;
  while (sitev <= PSIMAX) {
    h = max(fabs(sitev)*HPSIOPSI,HPSIMIN);

    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();

  sitepack.push_back(xsite);
  xsite.clear();
  // ----------------------------------

  vector<double> params = {MAXSTEP,TOL,2,RHOC,(double)sitepack[0].size(),TIMESTEP,NMAX,DELTAN,RECURSION}; // set parameters

  vector< vector<double> > xpi = {{PHIIN,PSIIN}}; // set i.c. for inflationary trajectories
  
  StocDeltaN sdn(MODEL,sitepack,xpi,0,params); // declare the system
  
  //sdn.sample(); // obtain 1 sample path
  //sdn.sample_logplot(); // plot obtained sample path
  
  sdn.solve(); // solve PDE & SDE to obtain power spectrum
  sdn.f_logplot(0); // show plot of <N>
  sdn.f_logplot(1); // show plot of <delta N^2>
  sdn.calP_plot(); // show plot of power spectrum of zeta


  // ---------- stop stop watch ----------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // -------------------------------------
}


// ---------- Lagrangian params. X[0]=phi, X[1]=psi ----------

double StocDeltaN::V(vector<double> &X)
{
  return LAMBDA4*((1-X[1]*X[1]/MM/MM)*(1-X[1]*X[1]/MM/MM)
		  + 2*X[0]*X[0]*X[1]*X[1]/PHIC/PHIC/MM/MM + (X[0]-PHIC)/MU1
		  - (X[0]-PHIC)*(X[0]-PHIC)/MU2/MU2);
}

double StocDeltaN::VI(vector<double> &X, int I) // \partial_I V
{
  if (I == 0) {
    return LAMBDA4*(1./MU1 - 2*(X[0]-PHIC)/MU2/MU2 + 4*X[0]*X[1]*X[1]/PHIC/PHIC/MM/MM);
  } else {
    return LAMBDA4*(4*X[0]*X[0]*X[1]/PHIC/PHIC/MM/MM
		    - 4*X[1]*(1-X[1]*X[1]/MM/MM)/MM/MM);
  }
}

double StocDeltaN::metric(vector<double> &X, int I, int J) // G_IJ
{
  if (I == J) {
    return 1;
  } else {
    return 0;
  }
}

double StocDeltaN::inversemetric(vector<double> &X, int I, int J) // G^IJ
{
  return metric(X,I,J);
}

double StocDeltaN::affine(vector<double> &X, int I, int J, int K) // \Gamma^I_JK
{
  return 0;
}

double StocDeltaN::derGamma(vector<double> &X, int I, int J, int K, int L) // \partial_L \Gamma^I_JK
{
  return 0;
}
