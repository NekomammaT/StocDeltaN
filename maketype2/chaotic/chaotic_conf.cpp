#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "chaotic_conf"

#define XMIN 0
#define XMAX 20
#define HX 1e-2

#define MAXSTEP 100000
#define TOL 1e-10

#define MPHI (1e-5)

#define RHOC (MPHI*MPHI)

#define RECURSION 100
#define XIN 13
#define TIMESTEP (1e-2)

#define DELTAN 0.1
#define NMAX 40


int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6; // start stop watch


  double h = HX, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > xsite;
  vector< vector< vector<double> > > sitepack;
  while (sitev <= XMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  sitepack.push_back(xsite);

  vector<double> params = {MAXSTEP,TOL,2,RHOC,(double)sitepack[0].size(),TIMESTEP,NMAX,DELTAN,
			   RECURSION};

  vector< vector<double> > xpi = {{XIN}};

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
  return 1./2*MPHI*MPHI*X[0]*X[0];
}

double StocDeltaN::VI(vector<double> &X, int I)
{
  return MPHI*MPHI*X[0];
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
