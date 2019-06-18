#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "USR_log"

#define XMIN -13
#define XMAX -1
#define HX 0.1
#define PMIN -(1e-12)
#define PMAX 1e-22
#define HPMIN (1e-23)
#define HPOV (1./20)

#define MAXSTEP 10000
#define TOL 1e-10

#define V0 (1e-20)
#define PHIF 0.000048005

#define RHOC V0

#define RECURSION 100
#define PDIN (1e-11)
#define PIN (-PDIN*PDIN/3./sqrt(V0/3.))
#define TIMESTEP (1e-4)
#define PHIINF (PDIN/3./sqrt(V0/3.))
#define XIN log(PHIINF)

#define DELTAN 0.01
#define NMAX 3


int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;

  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;

  double h = HX, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > xpsite;
  vector< vector< vector<double> > > sitepack;
  while (sitev <= XMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xpsite.push_back(site);
  sitepack.push_back(xpsite);
  site.clear();
  xpsite.clear();

  sitev = PMIN;
  while (sitev <= PMAX) {
    h = max(fabs(sitev)*HPOV,HPMIN);

    site.push_back(sitev);
    sitev += h;
  }
  xpsite.push_back(site);
  sitepack.push_back(xpsite);
  site.clear();
  xpsite.clear();

  vector<double> params = {MAXSTEP,TOL,2,RHOC,(double)sitepack[0].size(),TIMESTEP,NMAX,DELTAN,RECURSION};
  
  vector< vector<double> > xpi = {{XIN},{PIN}};

  StocDeltaN sdn(MODEL,sitepack,xpi,0,params);
  
  //sdn.sample();
  //sdn.sample_logplot();
  
  sdn.solve();
  sdn.f_loglogplot(0);
  sdn.f_loglogplot(1);
  sdn.calP_plot();

  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
  cout << after - before << " sec." << endl;
}



double StocDeltaN::V(vector<double> &X)
{
  if (X[0] > log(PHIF)) {
    return V0;
  } else {
    return 0;
  }
}

double StocDeltaN::VI(vector<double> &X, int I)
{
  return 0;
}

double StocDeltaN::metric(vector<double> &X, int I, int J)
{
  return exp(2*X[0]);
}

double StocDeltaN::inversemetric(vector<double> &X, int I, int J)
{
  return exp(-2*X[0]);
}

double StocDeltaN::affine(vector<double> &X, int I, int J, int K)
{
  return 1;
}

double StocDeltaN::derGamma(vector<double> &X, int I, int J, int K, int L)
{
  return 0;
}
