#include "../source/StocDeltaN_conf.hpp"
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
  vector< vector<double> > sitepack;
  while (sitev <= XMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  sitepack.push_back(site);
  site.clear();
  
  vector<double> xi = {XIN};
  
  StocDeltaN sdn(MODEL,sitepack,RHOC,xi,0,MAXSTEP,TOL,RECURSION,
		 TIMESTEP,NMAX,DELTAN);
  
  //sdn.sample();
  
  sdn.solve();

  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}


// --------------------- sample ---------------------------

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

double StocDeltaN::Dphi(vector<double> &X, int I)
{
  double Dphi = 0;

  for (int J=0; J<dim; J++) {
    Dphi -= inversemetric(X,I,J)*VI(X,J)/V(X);
  }

  return Dphi;
}

double StocDeltaN::Dphiphi(vector<double> &X, int I, int J)
{
  return V(X)/12./M_PI/M_PI * inversemetric(X,I,J);
}

double StocDeltaN::PhiNoise(vector<double> &X, int I, int alpha)
{
  return sqrt(V(X)/3.)/2./M_PI * vielbein(X,I,alpha);
}

