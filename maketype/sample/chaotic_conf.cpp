#include "../source/StocDeltaN_conf.hpp"
#include <sys/time.h>

#define MODEL "chaotic_conf" // model name

// ---------- box size & step h ------------
#define XMIN 0
#define XMAX 20
#define HX 1e-2
// -----------------------------------------

// ---------- for PDE ----------------------
#define MAXSTEP 100000 // max recursion
#define TOL 1e-10 // tolerance
// -----------------------------------------

#define MPHI (1e-5) // potential parameter

#define RHOC (MPHI*MPHI) // end of inflation

// ---------- for SDE ----------------------
#define RECURSION 100 // recursion for power spectrum
#define XIN 13 // i.c. for phi
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
  vector< vector<double> > sitepack;
  while (sitev <= XMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  sitepack.push_back(site);
  site.clear();
  // -----------------------------------------
  
  vector<double> xi = {XIN}; // set i.c. for sample paths
  
  StocDeltaN sdn(MODEL,sitepack,RHOC,xi,0,MAXSTEP,TOL,RECURSION,
		 TIMESTEP,NMAX,DELTAN); // declare the system
  
  //sdn.sample(); // show 1 sample path
  //sdn.sample_plot(); // plot that sample path
  
  sdn.solve(); // solve PDE & SDE and obtain power spectrum
  sdn.f1_plot(); // show plot of <N>
  sdn.g2_plot(); // show plot of <delta N^2>
  sdn.calP_plot(); // show plot of power spectrum of zeta

  // --------- stop stopwatch ----------------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // -----------------------------------------
}


// ---------- Lagrangian params. and diff. coeff.  X[0]=phi ----------

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

double StocDeltaN::Dphi(vector<double> &X, int I) // D^I
{
  double Dphi = 0;

  for (int J=0; J<dim; J++) {
    Dphi -= inversemetric(X,I,J)*VI(X,J)/V(X);
  }

  return Dphi;
}

double StocDeltaN::Dphiphi(vector<double> &X, int I, int J) // D^IJ
{
  return V(X)/12./M_PI/M_PI * inversemetric(X,I,J);
}

double StocDeltaN::PhiNoise(vector<double> &X, int I, int alpha) // g^I_a
{
  return sqrt(V(X)/3.)/2./M_PI * vielbein(X,I,alpha);
}

