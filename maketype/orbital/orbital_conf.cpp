#include "../source/StocDeltaN_conf.hpp"
#include <sys/time.h>

#define MODEL "orbital_conf" // model name

// ---------- box size & step h ------------
#define RHOMIN 0.01
#define RHOMAX 4
#define THETAMIN 0
#define THETAMAX 15
#define HRHO (1e-2)
#define HTHETA (1e-1)
// -----------------------------------------

// ---------- for PDE ----------------------
#define MAXSTEP 100000 // max recursion
#define TOL 1e-10 // tolerance
// -----------------------------------------

// ---------- potential parameter ----------
#define MM (1e-5) 
// -----------------------------------------

#define RHOC (MM*MM) // end of inflation

// ---------- for SDE ----------------------
#define RECURSION 100 // recursion for power spectrum
#define RHOIN 2 // i.c. for phi
#define THETAIN 10 // i.c. for psi
#define TIMESTEP (1e-2) // time step: delta N
// -----------------------------------------

// ---------- for power spectrum -----------
#define DELTAN 0.1 // calc. PS every DELTAN e-folds
#define NMAX 60 // calc. PS for 0--NMAX e-folds
// -----------------------------------------



int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6; // start stop watch
  
  // ---------- set box step h ---------------
  double h = HRHO, sitev = RHOMIN;
  vector<double> site;
  vector< vector<double> > sitepack;
  while (sitev <= RHOMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  sitepack.push_back(site);
  site.clear();

  h = HTHETA, sitev = THETAMIN;
  while (sitev <= THETAMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  sitepack.push_back(site);
  site.clear();
  // -----------------------------------------
  
  vector<double> xi = {RHOIN,THETAIN}; // set i.c. for sample paths
  
  StocDeltaN sdn(MODEL,sitepack,RHOC,xi,0,MAXSTEP,TOL,RECURSION,
		 TIMESTEP,NMAX,DELTAN); // declare the system
  
  //sdn.sample(); // show 1 sample path
  //sdn.sample_plot();
  
  sdn.solve(); // solve PDE & SDE and obtain power spectrum
  sdn.f1_plot();
  sdn.g2_plot();
  sdn.calP_plot();

  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}


// ---------- Lagrangian params. and diff. coeff.  X[0]=rho, X[1]=theta ----------

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

double StocDeltaN::Dphi(vector<double> &X, int I) // D^I
{
  double Dphi = 0;

  for (int J=0; J<dim; J++) {
    Dphi -= inversemetric(X,I,J)*VI(X,J)/V(X);
    
    for (int K=0; K<dim; K++) {
      Dphi -= affine(X,I,J,K)*Dphiphi(X,J,K)/2.;
    }
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
