#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "hilltop_conf"

#define XMIN -21
#define XMAX 21
#define HMIN (1e-8)
#define HOV (1./20)

#define MAXSTEP 100000000
#define TOL 1e-10

#define LAMBDA 0.01
#define MU 20

#define RHOC (LAMBDA*LAMBDA*LAMBDA*LAMBDA*(sqrt(1+2*MU*MU)-1)/MU/MU)

#define RECURSION 100
#define XIN (1e-6)
#define TIMESTEP (1e-1)

#define DELTAN 10
#define NMAX 2800


int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6; // start stop watch

  double h, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > xsite;
  vector< vector< vector<double> > > sitepack;
  while (sitev <= XMAX) {
    h = max(fabs(sitev)*HOV,HMIN);
    
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  sitepack.push_back(xsite);

  vector<double> params = {MAXSTEP,TOL,2,RHOC,(double)sitepack[0].size(),
			   TIMESTEP,NMAX,DELTAN,RECURSION};

  vector< vector<double> > xpi = {{XIN}};

  StocDeltaN sdn(MODEL,sitepack,xpi,0,params);

  sdn.solve();
  sdn.f_loglinearplot(0);
  sdn.f_loglinearplot(1);
  sdn.calP_plot();
  

  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}


// ------------------- user decision -----------------------
// ---------------------------------------------------------

double StocDeltaN::V(vector<double> &X)
{
  return LAMBDA*LAMBDA*LAMBDA*LAMBDA*(1-X[0]*X[0]/MU/MU);
}

double StocDeltaN::VI(vector<double> &X, int I)
{
  return -2*X[0]*LAMBDA*LAMBDA*LAMBDA*LAMBDA/MU/MU;
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
