#include "../source/JacobiPDE.hpp"

class Test: virtual public JacobiPDE
{
public:
  Test(){}
  Test(vector< vector< vector<double> > > &Site, vector<double> &Params);
  virtual double DI(int xp, int I, vector< vector<double> > &psv);
  virtual double DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv);
  virtual double CC(int num, vector< vector<double> > &psv, int func);
  virtual void BoundaryCondition();
};

int main()
{
  vector< vector< vector<double> > > site;
  vector< vector<double> > xsite;
  xsite.push_back({0,1,2,3,4,5,6,7,8,9});
  site.push_back(xsite);
  site.push_back(xsite);

  vector<double> params = {10000,1e-10,1};

  Test test(site,params);
  
  test.PDE_solve(0);
  test.export_fg("Laplace2.dat");
}

Test::Test(vector< vector< vector<double> > > &Site, vector<double> &Params):
  JacobiPDE(Site,Params)
{
  BoundaryCondition();
}

// ---------------------------------------------------------
/* 
solve (DI(xp,I) \partial_xpI + 1./2 DIJ(xpI,xpJ) \partial_xpI \partial_xpJ) f = CC
func swithes f.
*/
double Test::DI(int xp, int I, vector< vector<double> > &psv)
{
  return 0;
}

double Test::DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv)
{
  if (xpI == xpJ && I == J) {
    return 1;
  } else {
    return 0;
  }
}

double Test::CC(int num, vector< vector<double> > &psv, int func)
{
  return 0;
}
// ---------------------------------------------------------

void Test::BoundaryCondition()
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int number=0; number<volume; number++) {
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
	
	ind_m[xp][I]--;
	hm[number][xp][I] = hI[xp][I][index[xp][I]-1];

	ind_p[xp][I]++;
	hp[number][xp][I] = hI[xp][I][index[xp][I]];

	num_m[number][xp][I] = Ind2No(ind_m);
	num_p[number][xp][I] = Ind2No(ind_p);

	
	for (int xptemp=0; xptemp<xpdim; xptemp++) {
	  for (int J=0; J<Idim; J++) {
	    if (xp!=xptemp || I!=J) {
	      ind_pp = index;
	      ind_pm = index;
	      ind_mm = index;

	      ind_mm[xp][I]--;
	      ind_pp[xp][I]++;
	      ind_pm[xp][I]++;
	      
	      ind_pm[xptemp][J]--;
	      ind_mm[xptemp][J]--;
	      ind_pp[xptemp][J]++;

	      num_pp[number][xp][I][xptemp][J] = Ind2No(ind_pp);
	      num_pm[number][xp][I][xptemp][J] = Ind2No(ind_pm);
	      num_mm[number][xp][I][xptemp][J] = Ind2No(ind_mm);
	    }
	  }
	}
      }
    }

    if (index[0][0] == 0) {
      ff[0][number] = siteNo[1][0]-1-index[1][0];
      Omega[number] = false;
    } else if (index[0][0] == siteNo[0][0]-1) {
      ff[0][number] = index[1][0];
      Omega[number] = false;
    } else if (index[1][0] == 0) {
      ff[0][number] = siteNo[0][0]-1-index[0][0];
      Omega[number] = false;
    } else if (index[1][0] == siteNo[1][0]-1) {
      ff[0][number] = index[0][0];
      Omega[number] = false;
    } else {
      ff[0][number] = rand()%10;
      Omega[number] = true;
    }
    
  }

  for (int num=0; num<volume; num++) {
    cout << ff[0][num] << ' ';
    if (num%10 == 9) {
      cout << endl;
    }
  }
  cout << endl;
  
  for (int num=0; num<volume; num++) {
    cout << Omega[num] << ' ';
    if (num%10 == 9) {
      cout << endl;
    }
  }
  cout << endl;
}

// ---------------------------------------------------------
// ---------------------------------------------------------

