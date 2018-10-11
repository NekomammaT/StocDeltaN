#include "StocDeltaN_conf.hpp"

StocDeltaN::StocDeltaN(string Model,
		       vector< vector<double> > &Site, double Rhoc,
		       vector<double> &Xi, double T0,
		       double Maxstep, double Tol, int Recursion,
		       double Timestep, double NNmax, double DeltaN):
  JacobiPDE(Site,Rhoc), SRKintegrater(Xi,T0)
{
  model = Model;
  maxstep = Maxstep;
  tol = Tol;
  recursion = Recursion;
  timestep = Timestep;
  Nmax = NNmax;
  deltaN = DeltaN;
}

void StocDeltaN::init_fn()
{
  for (int number=0; number<volume; number++) {
    for (int I=0; I<dim; I++) {
      FPoint[I] = No2X(number,I);
    }
    
    if (V(FPoint) < rhoc) {
      Omega[number] = false;
      f1[number] = 0;
      g2[number] = 0;
    } else {
      Omega[number] = true;
      f1[number] = rand()%10;
      g2[number] = (rand()%10)/10.;
    }
  }
}

void StocDeltaN::init_txp()
{
  t = t0;
  x = xi;
}

void StocDeltaN::solve()
{
  init_fn();
  
  cout << model << endl;
  cout << "total # of sites: " << volume << endl;
  
  PDE_solve(maxstep,tol,1); //solve f1
  PDE_solve(maxstep,tol,-2);
  
  string str = "Mn_" + model + ".dat";
  
  export_fg(str); // export f1 and g2 to file
  
  
  vector< vector<double> > dN2List[recursion];
  str = "traj_" + model + ".dat";
  ofstream trajfile(str);
  double dt = timestep;
  int recNo = 0;

  init_txp();

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<recursion; i++) {
    StocDeltaN ssdn = *this;
    vector< vector<double> > dN2data;
    double dN2, dataNo;
    
    while (ssdn.return_V() > rhoc) {
      ssdn.SRK2(dt);
      
      if (i == 0) {
	trajfile << ssdn.return_t() << ' ';
	for (int I=0; I<dim; I++) {
	  trajfile << ssdn.return_phi(I) << ' ';
	}
	trajfile << ssdn.return_intf1() << ' ' << ssdn.return_intg2() << endl;
      }
      
      dN2data.push_back({ssdn.return_intf1(),ssdn.return_intg2()});
    }
    
    for (double N=0; N<Nmax; N+=deltaN) {
      dN2 = 0;
      dataNo = 0;
      
      for (int list=0; list<dN2data.size(); list++) {
	if (N <= dN2data[list][0] && dN2data[list][0] < N+deltaN) {
	  dN2 += dN2data[list][1];
	  dataNo++;
	}
      }
      
      if (dataNo != 0) {
	dN2 /= dataNo;
      }
      dN2List[i].push_back({N,dN2});
    }

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      recNo++;
      cout << "\r" << recNo << "/" << recursion << flush;
    }
  }
  cout << endl;
  
  double meandN2, predN2 = 0;
  str = "calP_" + model + ".dat";
  ofstream calPfile(str);

  for (int list=0; list<dN2List[0].size(); list++) {
    meandN2 = 0;
    recNo = 0;

    for (int i=0; i<recursion; i++) {
      meandN2 += dN2List[i][list][1];
      if (dN2List[i][list][1] != 0) {
	recNo++;
      }
    }

    meandN2 /= recNo;
    calPfile << dN2List[0][list][0] << ' ' << meandN2 << ' ' << (meandN2-predN2)/deltaN << endl;
    predN2 = meandN2;
  }
}

void StocDeltaN::sample()
{
  init_txp();

  vector<double> xmax = x, xmin = x;
  double Vi = return_V();

  string str = "sample_" + model + ".dat";
  ofstream ofs(str);
  double dt = timestep;

  while (return_V() > rhoc){
    SRK2(dt);

    for (int I=0; I<dim; I++) {
      if (x[I] < xmin[I]) {
	xmin[I] = x[I];
      }
      if (x[I] > xmax[I]) {
	xmax[I] = x[I];
      }
    }

    ofs << setprecision(6) << return_t() << ' ';
    for (int I=0; I<dim; I++) {
      ofs << return_phi(I) << ' ';
    }
    ofs << setprecision(17)
	<< return_V() << endl;
  }

  cout << "[xi, xf, xmin, xmax]: " << endl;
  for (int I=0; I<dim; I++) {
    cout << "[" << xi[I] << ", " << x[I] << ", " << xmin[I] << ", " << xmax[I] << "] " << endl;
  }
  cout << endl;

  cout << "N = " << return_t() << endl;
  cout << setprecision(17)
       << "Vi = " << Vi << ",  Vf = " <<return_V() << endl;
}

double StocDeltaN::return_intf1()
{
  return Interpolation_f(x,f1);
}

double StocDeltaN::return_intg2()
{
  return Interpolation_f(x,g2);
}
