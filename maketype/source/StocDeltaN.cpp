#include "StocDeltaN.hpp"

StocDeltaN::StocDeltaN(string Model,
		       vector< vector<double> > Site[], double Rhoc,
		       vector<double> &Xi, vector<double> &Pi, double T0, int NoiseDim,
		       double Maxstep, double Tol, int Recursion,
		       double Timestep, double NNmax, double DeltaN):
  JacobiPDE(Site,Rhoc), SRKintegrater(Xi,Pi,T0,NoiseDim)
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
      FPoint[0][I] = No2X(number,I);
      FPoint[1][I] = No2P(number,I);
    }
    
    if (V(FPoint[0]) < rhoc) {
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
  p = pi;
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
    
    while (3*ssdn.return_H()*ssdn.return_H() > rhoc) {
      ssdn.SRK2(dt);
      
      if (i == 0) {
	trajfile << ssdn.return_t() << ' ';
	for (int I=0; I<dim; I++) {
	  trajfile << ssdn.return_phi(I) << ' ';
	}
	for (int I=0; I<dim; I++) {
	  trajfile << ssdn.return_pi(I) << ' ';
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

  vector<double> xmax = x, xmin = x, pmax = p, pmin = p;
  double Hi = return_H();

  string str = "sample_" + model + ".dat";
  ofstream ofs(str);
  double dt = timestep;

  while (3*return_H()*return_H() > rhoc){ //return_e1() < 1) {
    SRK2(dt);

    for (int I=0; I<dim; I++) {
      if (x[I] < xmin[I]) {
	xmin[I] = x[I];
      }
      if (x[I] > xmax[I]) {
	xmax[I] = x[I];
      }
      if (p[I] < pmin[I]) {
	pmin[I] = p[I];
      }
      if (p[I] > pmax[I]) {
	pmax[I] = p[I];
      }
    }

    ofs << setprecision(6) << return_t() << ' ';
    for (int I=0; I<dim; I++) {
      ofs << return_phi(I) << ' ';
    }
    for (int I=0; I<dim; I++) {
      ofs << return_pi(I) << ' ';
    }
    ofs << setprecision(17)
	<< return_H() << endl;
  }

  cout << "[xi, xf, xmin, xmax]: " << endl;
  for (int I=0; I<dim; I++) {
    cout << "[" << xi[I] << ", " << x[I] << ", " << xmin[I] << ", " << xmax[I] << "] " << endl;
  }
  cout << endl;

  cout << "[pi, pf, pmin, pmax]: " << endl;
  for (int I=0; I<dim; I++) {
    cout << "[" << pi[I] << ", " << p[I] << ", " << pmin[I] << ", " << pmax[I] << "] " << endl;
  }
  cout << endl;

  cout << "N = " << return_t() << endl;
  cout << setprecision(17)
       << "Hi = " << Hi << ",  Hf = " <<return_H() << endl;
  cout << "eH = " << setprecision(6) << return_e1() << endl;
}

double StocDeltaN::return_intf1()
{
  return Interpolation_f(x,p,f1);
}

double StocDeltaN::return_intg2()
{
  return Interpolation_f(x,p,g2);
}
