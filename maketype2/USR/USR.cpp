#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "USR"

#define XMIN 0.04
#define XMAX 0.2
#define PMIN -(1e-10)
#define PMAX 0
#define HXMIN (1e-8)
#define HXOV (1./20)
#define HPMIN (1e-20)
#define HPOV (1./20)

#define MAXSTEP 100
#define TOL 1e-10

#define V0 (1e-20)
#define XF 0.042312978

#define RHOC V0

#define RECURSION 100
#define XIN 0.1 
#define PIN -(1e-11)
#define TIMESTEP (1e-3)

#define DELTAN 0.01
#define NMAX 6


int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;

  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;

  double h, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > xpsite;
  vector< vector< vector<double> > > sitepack;
  while (sitev <= XMAX) {
    h = max(fabs(sitev-XF)*HXOV,HXMIN);

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
  //sdn.sample_loglogplot();
  
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
  if (X[0] > XF) {
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

double StocDeltaN::derGamma(vector<double> &X, int I, int J, int K, int L)
{
  return 0;
}
