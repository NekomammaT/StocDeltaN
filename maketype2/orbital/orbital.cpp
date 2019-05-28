#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "orbital"

#define RHOMIN 0.99
#define RHOMAX 1.01
#define HRHO (1e-3)
#define THETAMIN 1
#define THETAMAX 13
#define HTHETA 1
#define RPIMIN -(1e-10)
#define RPIMAX (1e-10)
#define HRPI (1e-11)
#define TPIMIN -(9e-6)
#define TPIMAX -(7e-6)
#define HTPI (1e-6)

#define MAXSTEP 100000
#define TOL 1e-10

#define MM (1e-5)

#define RHOC (MM*MM)

#define RECURSION 100
#define RHOIN 1
#define RPIIN 0
#define THETAIN 10
#define TPIIN -sqrt(2./3)*MM
#define TIMESTEP (1e-2)

#define DELTAN 0.1
#define NMAX 24


int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6; // start stop watch


  double h = HRHO, sitev = RHOMIN;
  vector<double> site;
  vector< vector<double> > xpsite;
  vector< vector< vector<double> > > sitepack;
  while (sitev <= RHOMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xpsite.push_back(site);
  site.clear();

  h = HTHETA, sitev = THETAMIN;
  while (sitev <= THETAMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xpsite.push_back(site);
  sitepack.push_back(xpsite);
  site.clear();
  xpsite.clear();

  h = HRPI, sitev = RPIMIN;
  while (sitev <= RPIMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xpsite.push_back(site);
  site.clear();

  h = HTPI, sitev = TPIMIN;
  while (sitev <= TPIMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xpsite.push_back(site);
  sitepack.push_back(xpsite);
  site.clear();
  xpsite.clear();

  vector<double> params = {MAXSTEP,TOL,2,RHOC,(double)sitepack[0].size(),TIMESTEP,NMAX,DELTAN,
			   RECURSION};

  vector< vector<double> > xpi = {{RHOIN,THETAIN},{RPIIN,TPIIN}};

  StocDeltaN sdn(MODEL,sitepack,xpi,0,params);
  
  //sdn.sample();
  //sdn.sample_plot();
  
  sdn.solve();
  sdn.calP_plot();

  
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}


// ------------------- user decision -----------------------
// ---------------------------------------------------------

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

double StocDeltaN::derGamma(vector<double> &X, int I, int J, int K, int L)
{
  if (L == 0) {
    if (I == 0 && J == 1 && K == 1) {
      return -1;
    } else if (I == 1 && J != K) {
      return -1./X[0]/X[0];
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

