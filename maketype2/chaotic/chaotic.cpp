#include "../source/JacobiPDE.hpp"
#include "../source/SRK32.hpp"
#include <sys/time.h>


class Chaotic: virtual public JacobiPDE, virtual public SRKintegrater
{
protected:
  int xpdim,Idim;
  
public:
  Chaotic(){}
  Chaotic(vector< vector< vector<double> > > &Site, vector<double> &Params,
	  vector< vector<double> > &XPi, double T0, int NoiseDim);
  virtual double V(vector<double> &X);
  virtual double VI(vector<double> &X, int I);
  virtual double metric(vector<double> &X, int I, int J);
  virtual double inversemetric(vector<double> &X, int I, int J);
  virtual double affine(vector<double> &X, int I, int J, int K);
  virtual double derGamma(vector<double> &X, int I, int J, int K, int L); // Gamma^I_{JK,L}
};


#define XMIN 0
#define XMAX 15
#define PMIN (-1e-4)
#define PMAX 0
#define HX 0.01
#define HPMIN (1e-6)
#define HPOV (1./20)

#define MAXSTEP 10000
#define TOL (1e-10)

#define MPHI (1e-5)

#define RHOC (MPHI*MPHI)

#define PHIIN 13
#define PIN -(5e-6)
#define TIMESTEP (1e-2)


int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6; // start stop watch


  double h = HX, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > xpsite;
  vector< vector< vector<double> > > sitepack;
  while(sitev <= XMAX) {
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

  vector<double> params = {MAXSTEP,TOL,2,RHOC};

  vector< vector<double> > xpi = {{PHIIN},{PIN}};

  Chaotic chaotic(sitepack,params,xpi,0,xpi[0].size());
  ofstream trajfile("traj_chaotic.dat");

  chaotic.PDE_solve(0);
  chaotic.PDE_solve(1);
  chaotic.export_fg("Mn_chaotic.dat");

  while (3*chaotic.return_H()*chaotic.return_H() >= RHOC) {
    chaotic.SRK2(TIMESTEP);
    
    trajfile << chaotic.return_t() << ' '
	     << chaotic.return_xp(0,0) << ' '
	     << chaotic.return_xp(1,0) << endl;
  }


  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}



Chaotic::Chaotic(vector< vector< vector<double> > > &Site, vector<double> &Params,
		 vector< vector<double> > &XPi, double T0, int NoiseDim):
  JacobiPDE(Site,Params), SRKintegrater(XPi,T0,NoiseDim)
{
  xpdim = JacobiPDE::xpdim;
  Idim = JacobiPDE::Idim;
  BoundaryCondition();
}

// ------------------- user decision -----------------------
// ---------------------------------------------------------

double Chaotic::V(vector<double> &X)
{
  return 1./2*MPHI*MPHI*X[0]*X[0];
}

double Chaotic::VI(vector<double> &X, int I)
{
  return MPHI*MPHI*X[0];
}

double Chaotic::metric(vector<double> &X, int I, int J)
{
  return 1;
}

double Chaotic::inversemetric(vector<double> &X, int I, int J)
{
  return 1;
}

double Chaotic::affine(vector<double> &X, int I, int J, int K)
{
  return 0;
}

double Chaotic::derGamma(vector<double> &X, int I, int J, int K, int L)
{
  return 0;
}

// ---------------------------------------------------------
// ---------------------------------------------------------

