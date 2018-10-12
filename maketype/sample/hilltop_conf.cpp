#include "../source/StocDeltaN_conf.hpp"
#include <sys/time.h>

#define MODEL "hilltop_conf" // model name

// ---------- box size & step h ------------
#define XMIN -21
#define XMAX 21
#define HMIN (1e-8)
#define HOV (1./20) // h/|phi|
// -----------------------------------------

// ---------- for PDE ----------------------
#define MAXSTEP 100000000 // max recursion
#define TOL 1e-10 // tolerance
// -----------------------------------------

// ---------- potential parameter ----------
#define LAMBDA 0.01
#define MU 20
// -----------------------------------------

#define RHOC (LAMBDA*LAMBDA*LAMBDA*LAMBDA*(sqrt(1+2*MU*MU)-1)/MU/MU) // end of inflation

// ---------- for SDE ----------------------
#define RECURSION 100 // recursion for power spectrum
#define XIN (1e-6) // i.c. for phi
#define TIMESTEP (1e-1) // time step: delta N
// -----------------------------------------

// ---------- for power spectrum -----------
#define DELTAN 10 // calc. PS every DELTAN e-folds
#define NMAX 2800 // calc. PS for 0--NMAX e-folds
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
  double h, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > sitepack;
  while (sitev <= XMAX) {
    h = max(fabs(sitev)*HOV,HMIN);

    site.push_back(sitev);
    sitev += h;
  }
  sitepack.push_back(site);
  site.clear();
  // ------------------------------------------

  vector<double> xi = {XIN}; // set i.c. for sample paths
  
  StocDeltaN sdn(MODEL,sitepack,RHOC,xi,0,MAXSTEP,TOL,RECURSION,
		 TIMESTEP,NMAX,DELTAN); // declare the system
  
  //sdn.sample(); // show 1 sample path
  
  sdn.solve(); // solve PDE & SDE and obtain power spectrum

  // --------- stop stopwatch ----------------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // -----------------------------------------
}


// ---------- Lagrangian params. and diff. coeff.  X[0]=phi ----------

double StocDeltaN::V(vector<double> &X)
{
  return LAMBDA*LAMBDA*LAMBDA*LAMBDA*(1-X[0]*X[0]/MU/MU);
}

double StocDeltaN::VI(vector<double> &X, int I) // \partial_I V
{
  return -2*X[0]*LAMBDA*LAMBDA*LAMBDA*LAMBDA/MU/MU;
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

