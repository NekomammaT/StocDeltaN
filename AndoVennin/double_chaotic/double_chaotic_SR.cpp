#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "double_chaotic_SR" // model name

// ---------- box size & step h --------------
#define PHIMIN -5
#define PHIMAX 20
#define PSIMIN 0
#define PSIMAX 20
#define HPHI 0.1
#define HPSI 0.1
// -------------------------------------------

// ---------- for PDE ------------------------
#define MAXSTEP 100000 // max recursion
#define TOL 1e-10 // tolerance
// -------------------------------------------

// ---------- potential parameter ------------
#define MPHI (1e-5)
#define MPSI (MPHI/9.)
// -------------------------------------------

#define RHOC (MPSI*MPSI) // energy density at the end of inflation

// ---------- for SDE ------------------------
#define RECURSION 100 // recursion for power spectrum
#define PHIIN 13 // i.c. for phi
#define PSIIN 13 // i.c. for psi
#define TIMESTEP (1e-2) // time step : delta N
// -------------------------------------------

// ---------- for power spectrum -------------
#define DELTAN 1 // calc. PS every DELTAN e-folds
#define NMAX 80 // calc. PS for 0--NMAX e-folds
// -------------------------------------------

//#define NCUT 300


int main(int argc, char** argv)
{
  // ---------------- start stop watch --------------
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  // ------------------------------------------------


  double h = HPHI, sitev = PHIMIN;
  vector<double> site;
  vector< vector<double> > xsite;
  vector< vector< vector<double> > > sitepack;
  while(sitev <= PHIMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();

  h = HPSI, sitev = PSIMIN;
  while(sitev <= PSIMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();

  sitepack.push_back(xsite);
  xsite.clear();

  vector<double> params = {MAXSTEP,TOL,2,RHOC,(double)sitepack[0].size(),TIMESTEP,
    NMAX,DELTAN,RECURSION //,NCUT
  };

  vector< vector<double> > xpi = {{PHIIN,PSIIN}};

  StocDeltaN sdn(MODEL,sitepack,xpi,0,params);

  //sdn.sample(); 
  //sdn.sample_plot();

  sdn.solve();
  sdn.f_plot(0);
  sdn.f_plot(1);
  sdn.calP_plot();

  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}


// ------------------- user decision -----------------------
// ---------------------------------------------------------

double StocDeltaN::V(vector<double> &X)
{
  return 1./2*MPHI*MPHI*X[0]*X[0] + 1./2*MPSI*MPSI*X[1]*X[1];
}

double StocDeltaN::VI(vector<double> &X, int I)
{
  if (I == 0) {
    return MPHI*MPHI*X[0];
  } else {
    return MPSI*MPSI*X[1];
  }
}

double StocDeltaN::metric(vector<double> &X, int I, int J)
{
  if (I == J) {
    return 1;
  } else {
    return 0;
  }
}

double StocDeltaN::inversemetric(vector<double> &X, int I, int J)
{
  return metric(X,I,J);
}

double StocDeltaN::affine(vector<double> &X, int I, int J, int K)
{
  return 0;
}

double StocDeltaN::derGamma(vector<double> &X, int I, int J, int K, int L)
{
  return 0;
}
