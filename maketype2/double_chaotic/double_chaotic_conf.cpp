#include "../source/JacobiPDE.hpp"
#include <sys/time.h>

#define PHIMIN -5
#define PHIMAX 20
#define PSIMIN 0
#define PSIMAX 20
#define HPHI 0.1
#define HPSI 0.1

#define MAXSTEP 100000
#define TOL 1e-10

#define MPHI (1e-5)
#define MPSI (MPHI/9.)

#define RHOC (MPSI*MPSI)


class DoubleChaotic: virtual public JacobiPDE
{
public:
  DoubleChaotic(){}
  DoubleChaotic(vector< vector< vector<double> > > &Site, vector<double> &Params);
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


  double h = HPHI, sitev = PHIMIN;
  vector<double> site;
  vector< vector<double> > xsite;
  vector< vector< vector<double> > > sitepack;
  while(sitev <= PHIMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();

  h = HPSI, sitev = PSIMIN;
  while(sitev <= PSIMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();

  sitepack.push_back(xsite);
  xsite.clear();

  vector<double> params = {MAXSTEP,TOL,2,RHOC};

  DoubleChaotic dc(sitepack,params);

  dc.PDE_solve(0);
  dc.PDE_solve(1);
  dc.export_fg("double_chaotic_conf.dat");


  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}


DoubleChaotic::DoubleChaotic(vector< vector< vector<double> > > &Site, vector<double> &Params):
  JacobiPDE(Site,Params)
{
  BoundaryCondition();
}

double DoubleChaotic::V(vector<double> &X)
{
  return 1./2*MPHI*MPHI*X[0]*X[0] + 1./2*MPSI*MPSI*X[1]*X[1];
}

double DoubleChaotic::VI(vector<double> &X, int I)
{
  if (I == 0) {
    return MPHI*MPHI*X[0];
  } else {
    return MPSI*MPSI*X[1];
  }
}

double DoubleChaotic::metric(vector<double> &X, int I, int J)
{
  if (I == J) {
    return 1;
  } else {
    return 0;
  }
}

double DoubleChaotic::inversemetric(vector<double> &X, int I, int J)
{
  return metric(X,I,J);
}

double DoubleChaotic::affine(vector<double> &X, int I, int J, int K)
{
  return 0;
}

double DoubleChaotic::DI(int xp, int I, vector< vector<double> > &psv, int func)
{
  double DI = 0;

  for (int J=0; J<Idim; J++) {
    DI -= inversemetric(psv[0],I,J)*VI(psv[0],J)/V(psv[0]);
  }

  return DI;
}

double DoubleChaotic::DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv, int func)
{
  return V(psv[0])/12./M_PI/M_PI * inversemetric(psv[0],I,J);
}

void DoubleChaotic::BoundaryCondition()
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

