#include "../source/StocDeltaN_conf.hpp"
#include <sys/time.h>

#define MODEL "double_chaotic_conf" // model name

// ---------- box size & step h ------------
#define PHIMIN -5
#define PHIMAX 20
#define PSIMIN 0
#define PSIMAX 20
#define HPHI (1e-1)
#define HPSI (1e-1)
// -----------------------------------------

// ---------- for PDE ----------------------
#define MAXSTEP 100000 // max recursion
#define TOL 1e-10 // tolerance
// -----------------------------------------

// ---------- potential parameter ----------
#define MPHI (1e-5) 
#define MPSI (MPHI/9.)
// -----------------------------------------

#define RHOC (MPSI*MPSI) // end of inflation

// ---------- for SDE ----------------------
#define RECURSION 100 // recursion for power spectrum
#define PHIIN 13 // i.c. for phi
#define PSIIN 13 // i.c. for psi
#define TIMESTEP (1e-2) // time step: delta N
// -----------------------------------------

// ---------- for power spectrum -----------
#define DELTAN 0.1 // calc. PS every DELTAN e-folds
#define NMAX 80 // calc. PS for 0--NMAX e-folds
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
  double h = HPHI, sitev = PHIMIN;
  vector<double> site;
  vector< vector<double> > sitepack;
  while (sitev <= PHIMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  sitepack.push_back(site);
  site.clear();

  h = HPSI, sitev = PSIMIN;
  while (sitev <= PSIMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  sitepack.push_back(site);
  site.clear();
  // -----------------------------------------
  
  vector<double> xi = {PHIIN,PSIIN}; // set i.c. for sample paths
  
  StocDeltaN sdn(MODEL,sitepack,RHOC,xi,0,MAXSTEP,TOL,RECURSION,
		 TIMESTEP,NMAX,DELTAN); // declare the system
  
  sdn.sample(); // obtain 1 sample path
  //sdn.sample_plot(); // plot that sample path
  
  //sdn.solve(); // solve PDE & SDE and obtain power spectrum
  //sdn.f1_plot(); // show plot of <N>
  //sdn.g2_plot(); // show plot of <delta N^2>
  //sdn.calP_plot(); // show plot of power spectrum of zeta

  // --------- stop stopwatch ----------------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // -----------------------------------------
}


// ---------- Lagrangian params. and diff. coeff.  X[0]=phi, X[1]=psi ----------

double StocDeltaN::V(vector<double> &X)
{
  return 1./2*MPHI*MPHI*X[0]*X[0] + 1./2*MPSI*MPSI*X[1]*X[1];
}

double StocDeltaN::VI(vector<double> &X, int I) // \partial_I V
{
  if (I == 0) {
    return MPHI*MPHI*X[0];
  } else {
    return MPSI*MPSI*X[1];
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
