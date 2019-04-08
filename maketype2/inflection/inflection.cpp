#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "inflection"

#define XMIN 0.3
#define XMAX 6
#define PMIN -5e-5
#define PMAX 0
#define HXOV 0.01
#define HXMIN (1e-2)
#define HXMAX 0.1
#define HPOV (1./20)
#define HPMIN (1e-9)

#define MAXSTEP 10000
#define TOL 1e-10

#define LAMBDA 1
#define XI 2.3
#define PHIC 0.66
#define ALPHA (6*LAMBDA*PHIC/(3+XI*XI*PHIC*PHIC*PHIC*PHIC) - 4.3e-5)
#define M2 (LAMBDA*PHIC*PHIC*(3+XI*PHIC*PHIC)/(3+XI*XI*PHIC*PHIC*PHIC*PHIC))
#define V0 (1e-7)

#define HF (1.9e-5)
#define RHOC (3*HF*HF)

#define RECURSION 100
#define XIN 5
#define PIN -6e-6
#define TIMESTEP (1e-2)

#define DELTAN 0.01
#define NMAX 55


int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6; // start stop watch


  double h, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > xpsite;
  vector< vector< vector<double> > > sitepack;
  while (sitev <= XMAX) {
    h = min(max(fabs(sitev-PHIC)*HXOV, HXMIN), HXMAX);

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

  vector<double> params = {MAXSTEP,TOL,2,RHOC,(double)sitepack[0].size(),
			   TIMESTEP,NMAX,DELTAN,RECURSION};

  vector< vector<double> > xpi = {{XIN},{PIN}};

  StocDeltaN sdn(MODEL,sitepack,xpi,0,params);

  //sdn.sample();
  //sdn.sample_logplot();

  sdn.solve();
  sdn.f_logplot(0);
  sdn.f_logplot(1);
  sdn.calP_plot();


  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}


// ------------------- user decision -----------------------
// ---------------------------------------------------------

double StocDeltaN::V(vector<double> &X)
{
  return V0/12.*(6*M2*X[0]*X[0] - 4*ALPHA*X[0]*X[0]*X[0]
		 + 3*LAMBDA*X[0]*X[0]*X[0]*X[0])
    /(1+XI*X[0]*X[0])/(1+XI*X[0]*X[0]);
}

double StocDeltaN::VI(vector<double> &X, int I)
{
  return V0*X[0]*(M2*(3-3*XI*X[0]*X[0])
		  + X[0]*(3*LAMBDA*X[0]+ALPHA*(-3+XI*X[0]*X[0])))
    /3./(1+XI*X[0]*X[0])/(1+XI*X[0]*X[0])/(1+XI*X[0]*X[0]);
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
