#include "../source/StocDeltaN_conf.hpp"

#define MODEL "inflection_conf"

#define XMIN 0.3
#define XMAX 6
#define FINEXMIN 0.6
#define FINEXMAX 1
#define HX (1e-2)
#define FINEHX (1e-4)

#define MAXSTEP 100000
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
#define TIMESTEP (1e-2)

#define DELTAN 0.1
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
  vector< vector<double> > sitepack;
  while (sitev <= XMAX) {
    if (FINEXMIN <= sitev && sitev <= FINEXMAX) {
      h = FINEHX;
    } else {
      h = HX;
    }
    
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
  return V0/12.*(6*M2*X[0]*X[0] - 4*ALPHA*X[0]*X[0]*X[0] + 3*LAMBDA*X[0]*X[0]*X[0]*X[0])
    /(1+XI*X[0]*X[0])/(1+XI*X[0]*X[0]);
}

double StocDeltaN::VI(vector<double> &X, int I)
{
  return V0*X[0]*(M2*(3-3*XI*X[0]*X[0]) + X[0]*(3*LAMBDA*X[0]+ALPHA*(-3+XI*X[0]*X[0])))
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

