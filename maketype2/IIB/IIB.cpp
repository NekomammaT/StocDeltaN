#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "IIB"

#define XMIN 2
#define XMAX 12
#define FINEXMIN 3.55
#define FINEXMAX 3.8
#define PMIN -1e-5
#define PMAX 0
#define HXOV (1./20)
#define HXMAX 0.1
#define HXMIN (1e-4)
#define HPOV (1./20)
#define HPMIN (1e-10)

#define MAXSTEP 2
#define TOL 1e-10

#define AW 0.02
#define BW 1
#define CW 0.04
#define DW 0
#define GW (3.076278e-2)
#define RW (7.071067e-1)
#define VV 1000
#define W0 12.35
#define CUP 0.0382

#define HF (2e-6)
#define RHOC (3*HF*HF)

#define RECURSION 100
#define XIN 10
#define PIN -1e-10
#define TIMESTEP (1e-3)

#define DELTAN 0.01
#define NMAX 65


int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6; // start stop watch


  double h = HXMAX, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > xpsite;
  vector< vector< vector<double> > > sitepack;
  while (sitev <= XMAX) {
    if (FINEXMIN <= sitev && sitev <= FINEXMAX) {
      h = HXMIN;
    } else {
      h = min(max(min(fabs(sitev-FINEXMAX),fabs(sitev-FINEXMIN))*HXOV, HXMIN), HXMAX);
    }

    //cout << sitev << ' ';
    site.push_back(sitev);
    sitev += h;
  }
  //cout << endl;
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
  //sdn.f_logplot(0);
  //sdn.f_logplot(1);
  //sdn.calP_plot();


  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}


// ------------------- user decision -----------------------
// ---------------------------------------------------------

double StocDeltaN::V(vector<double> &X)
{
  return W0*W0/VV/VV/VV * (CUP/pow(VV,1./3) + AW/(exp(X[0]/sqrt(3))-BW)
			   - CW/exp(X[0]/sqrt(3))
			   + exp(2*X[0]/sqrt(3))/VV
			   *(DW-GW/(RW*exp(sqrt(3)*X[0])/VV + 1)) );
}

double StocDeltaN::VI(vector<double> &X, int I)
{
  return ((CW - (AW*exp((2*X[0])/sqrt(3)))/pow(BW - exp(X[0]/sqrt(3)),2) + 
	   exp(sqrt(3)*X[0])*((2*DW)/VV + 
			      (GW*(exp(sqrt(3)*X[0])*RW - 2*VV))/pow(exp(sqrt(3)*X[0])*RW + VV,2)))*
	  pow(W0,2))/(sqrt(3)*exp(X[0]/sqrt(3))*pow(VV,3));
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
