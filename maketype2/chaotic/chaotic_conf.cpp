#include "../source/JacobiPDE.hpp"
#include <sys/time.h>

#define XMIN 0
#define XMAX 20
#define HX 1e-2

#define MAXSTEP 100000
#define TOL 1e-10

#define MPHI (1e-5)

#define RHOC (MPHI*MPHI)


class Chaotic: virtual public JacobiPDE
{
public:
  Chaotic(){}
  Chaotic(vector< vector< vector<double> > > &Site, vector<double> &Params);
  virtual double V(vector<double> &X);
  virtual double VI(vector<double> &X, int I);
  virtual double metric(vector<double> &X, int I, int J);
  virtual double inversemetric(vector<double> &X, int I, int J);
  virtual double affine(vector<double> &X, int I, int J, int K);
  virtual double DI(int xp, int I, vector< vector<double> > &psv, int func);
  virtual double DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv, int func);
  virtual void BoundaryCondition();
};


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

  vector<double> params = {MAXSTEP,TOL,2,RHOC};
  
  Chaotic chaotic(sitepack,params);

  chaotic.PDE_solve(0);
  chaotic.PDE_solve(1);
  chaotic.export_fg("chaotic_conf.dat");
  

  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}


Chaotic::Chaotic(vector< vector< vector<double> > > &Site, vector<double> &Params):
  JacobiPDE(Site,Params)
{
  BoundaryCondition();
}

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

double Chaotic::DI(int xp, int I, vector< vector<double> > &psv, int func)
{
  double DI = 0;

  for (int J=0; J<Idim; J++) {
    DI -= inversemetric(psv[0],I,J)*VI(psv[0],J)/V(psv[0]);
  }

  return DI;
}

double Chaotic::DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv, int func)
{
  return V(psv[0])/12./M_PI/M_PI * inversemetric(psv[0],I,J);
}

void Chaotic::BoundaryCondition()
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int number=0; number<volume; number++) {
    vector< vector<double> > PSV0(xpdim, vector<double>(Idim,0));

    for (int xp=0; xp<xpdim; xp++) {
      for (int I=0; I<Idim; I++) {
	PSV0[xp][I] = No2PSV(number,xp,I);
      }
    }

    if (V(PSV0[0]) < rhoc) {
      Omega[number] = false;
      for (int func=0; func<funcNo; func++) {
	ff[func][number] = 0;
      }
    } else {
      Omega[number] = true;
      for (int func=0; func<funcNo; func++) {
	ff[func][number] = rand()%10;
      }
    }
    
    
    vector< vector<int> > index(xpdim, vector<int>(Idim,0));
    vector< vector<int> > ind_p, ind_m, ind_pp, ind_pm, ind_mm;

    for (int xp=0; xp<xpdim; xp++) {
      for (int I=0; I<Idim; I++) {
	index[xp][I] = No2Ind(number,xp,I);
      }
    }

    for (int xp=0; xp<xpdim; xp++) {
      for (int I=0; I<Idim; I++) {
	ind_p = index;
	ind_m = index;
	
	if (index[xp][I] == 0) {
	  ind_m[xp][I]++;
	  hm[number][xp][I] = hI[xp][I][index[xp][I]];
	} else {
	  ind_m[xp][I]--;
	  hm[number][xp][I] = hI[xp][I][index[xp][I]-1];
	}

	if (index[xp][I] == siteNo[xp][I]-1) {
	  ind_p[xp][I]--;
	  hp[number][xp][I] = hI[xp][I][index[xp][I]-1];
	} else {
	  ind_p[xp][I]++;
	  hp[number][xp][I] = hI[xp][I][index[xp][I]];
	}

	num_m[number][xp][I] = Ind2No(ind_m);
	num_p[number][xp][I] = Ind2No(ind_p);

	
	for (int xptemp=0; xptemp<xpdim; xptemp++) {
	  for (int J=0; J<Idim; J++) {
	    if (xp!=xptemp || I!=J) {
	      ind_pp = index;
	      ind_pm = index;
	      ind_mm = index;
	      
	      if (index[xp][I] == 0) {
		ind_mm[xp][I]++;
	      } else {
		ind_mm[xp][I]--;
	      }

	      if (index[xp][I] == siteNo[xp][I]-1) {
		ind_pp[xp][I]--;
		ind_pm[xp][I]--;
	      } else {
		ind_pp[xp][I]++;
		ind_pm[xp][I]++;
	      }

	      if (index[xptemp][J] == 0) {
		ind_pm[xptemp][J]++;
		ind_mm[xptemp][J]++;
	      } else {
		ind_pm[xptemp][J]--;
		ind_mm[xptemp][J]--;
	      }

	      if (index[xptemp][J] == siteNo[xptemp][J]-1) {
		ind_pp[xptemp][J]--;
	      } else {
		ind_pp[xptemp][J]++;
	      }

	      num_pp[number][xp][I][xptemp][J] = Ind2No(ind_pp);
	      num_pm[number][xp][I][xptemp][J] = Ind2No(ind_pm);
	      num_mm[number][xp][I][xptemp][J] = Ind2No(ind_mm);
	    }
	  }
	}
      }
    }
  }
}
