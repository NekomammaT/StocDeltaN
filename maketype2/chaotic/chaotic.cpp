#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "chaotic"

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

#define RECURSION 100
#define PHIIN 13
#define PIN -(5e-6)
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

  vector<double> params = {MAXSTEP,TOL,2,RHOC,(double)sitepack[0].size(),TIMESTEP,NMAX,DELTAN,
			   RECURSION};

  vector< vector<double> > xpi = {{PHIIN},{PIN}};

  StocDeltaN sdn(MODEL,sitepack,xpi,0,params);

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

double StocDeltaN::H(vector<double> &X, vector<double> &P)
{
  double rho = V(X);

  for (int I=0; I<X.size(); I++) {
    for (int J=0; J<X.size(); J++) {
      rho += 1./2*inversemetric(X,I,J)*P[I]*P[J];
    }
  }

  return sqrt(rho/3.);
}

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

// ---------------------------------------------------------
/* 
solve (DI(xp,I) \partial_xpI + 1./2 DIJ(xpI,xpJ) \partial_xpI \partial_xpJ) f = CC
func swithes f.
*/
double StocDeltaN::DI(int xp, int I, vector< vector<double> > &psv)
{
  double DI = 0;

  if (xpdim == 1) {
    for (int J=0; J<Idim; J++) {
      DI -= inversemetric(psv[0],I,J)*VI(psv[0],J)/V(psv[0]);
    }
  } else if (xpdim == 2) {
    if (xp == 0) {
      DI = 0;
      
      for (int J=0; J<Idim; J++) {
	DI += inversemetric(psv[0],I,J)*psv[1][J]/H(psv[0],psv[1]);
	
	for (int K=0; K<Idim; K++) {
	  DI -= 1./2*affine(psv[0],I,J,K)*DIJ(0,J,0,K,psv);
	}
      }
    } else {
      double Hubble = H(psv[0],psv[1]);
      DI = -3*psv[1][I]-VI(psv[0],I)/Hubble;
      
      for (int J=0; J<Idim; J++) {
	for (int K=0; K<Idim; K++) {
	  DI += affine(psv[0],J,I,K)*DIJ(0,K,1,J,psv);
	  
	  for (int S=0; S<Idim; S++) {
	    DI += affine(psv[0],S,I,J)*inversemetric(psv[0],J,K)*psv[1][K]*psv[1][S]/Hubble
	      +1./2*derGamma(psv[0],S,I,J,K)*psv[1][S]*DIJ(0,J,0,K,psv);
	    
	    for (int R=0; R<Idim; R++) {
	      DI -= 1./2*(affine(psv[0],R,J,K)*affine(psv[0],S,I,R)
			  + affine(psv[0],R,I,J)*affine(psv[0],S,K,R))
		*psv[1][S]*DIJ(0,J,0,K,psv);
	    }
	  }
	}
      }
    }
  }
  
  return DI;
}

double StocDeltaN::DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv)
{
  double DDIJ;
  
  if (xpdim == 1) {
    DDIJ = V(psv[0])/12./M_PI/M_PI * inversemetric(psv[0],I,J);
  } else if (xpdim == 2) {
    if (xpI == 0 && xpJ == 0) {
      DDIJ = H(psv[0],psv[1])*H(psv[0],psv[1])/4./M_PI/M_PI * inversemetric(psv[0],I,J);
    } else if (xpI == 1 && xpJ == 1) {
      DDIJ = 0;
      for (int K=0; K<Idim; K++) {
	for (int L=0; L<Idim; L++) {
	  for (int M=0; M<Idim; M++) {
	    for (int N=0; N<Idim; N++) {
	      DDIJ += affine(psv[0],K,I,L)*psv[1][K]*affine(psv[0],M,J,N)*psv[1][M]
		*DIJ(0,L,0,N,psv);
	    }
	  }
	}
      }
    } else if (xpI == 0) {
      DDIJ = 0;
      
      for (int K=0; K<Idim; K++) {
	for (int L=0; L<Idim; L++) {
	  DDIJ += affine(psv[0],K,J,L)*psv[1][K]*DIJ(0,I,0,L,psv);
	}
      }
    } else {
      DDIJ = DIJ(xpJ,J,xpI,I,psv);
    }
  }
  
  return DDIJ;
}

double StocDeltaN::CC(int num, vector< vector<double> > &psv, int func)
{
  double CC = 0;
  
  if (func == 0) {
    CC = -1;
  } else if (func == 1) {
    for (int xpI=0; xpI<xpdim; xpI++) {
      for (int I=0; I<Idim; I++) {
	for (int xpJ=0; xpJ<xpdim; xpJ++) {
	  for (int J=0; J<Idim; J++) {
	    CC -= DIJ(xpI,I,xpJ,J,psv)
	      *(ff[0][num_p[num][xpI][I]] - ff[0][num_m[num][xpI][I]])
	      *(ff[0][num_p[num][xpJ][J]] - ff[0][num_m[num][xpJ][J]])
	      /(hp[num][xpI][I]+hm[num][xpI][I])/(hp[num][xpJ][J]+hm[num][xpJ][J]);
	  }
	}
      }
    }
  }

  return CC;
}
// ---------------------------------------------------------

double StocDeltaN::gIa(int xp, int I, int alpha, vector< vector<double> > &psv)
{
  double ggIa;
  
  if (xpdim == 1) {
    ggIa = sqrt(V(psv[0])/3.)/2./M_PI * vielbein(psv,I,alpha);
  } else if (xpdim == 2) {
    if (xp == 0) {
      ggIa = H(psv[0],psv[1])/2./M_PI * vielbein(psv,I,alpha);
    } else if (xp == 1) {
      ggIa = 0;
      
      for (int J=0; J<Idim; J++) {
	for (int K=0; K<Idim; K++) {
	  ggIa += affine(psv[0],K,I,J)*psv[1][K]*gIa(0,J,alpha,psv);
	}
      }
    }
  }

  return ggIa;
}

void StocDeltaN::BoundaryCondition()
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

    if (xpdim == 1) {
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
    } else if (xpdim == 2) {
      if (3*H(PSV0[0],PSV0[1])*H(PSV0[0],PSV0[1]) < rhoc) {
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

// ---------------------------------------------------------
// ---------------------------------------------------------

