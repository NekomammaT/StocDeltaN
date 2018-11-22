#include "StocDeltaN_conf.hpp"
#include "matplotlibcpp.hpp"

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
  
  
  vector< vector< vector<double> > > dN2List(recursion);
  str = "traj_" + model + ".dat";
  ofstream trajfile(str);
  double dt = timestep;
  int recNo = 0;

  if (dim == 2) {
    x1traj.clear();
    x2traj.clear();
  }
  
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

	if (dim == 2) {
	  x1traj.push_back(ssdn.return_phi(0));
	  x2traj.push_back(ssdn.return_phi(1));
	}
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
    Ndata.push_back(dN2List[0][list][0]);
    calPdata.push_back((meandN2-predN2)/deltaN);
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

    if (dim == 1) {
      Ntraj.push_back(return_t());
      x1traj.push_back(return_phi(0));
    } else if (dim == 2) {
      x1traj.push_back(return_phi(0));
      x2traj.push_back(return_phi(1));
    }
  }

  cout << "[xi, xf, xmin, xmax]: " << endl;
  for (int I=0; I<dim; I++) {
    cout << "[" << xi[I] << ", " << x[I] << ", " << xmin[I] << ", " << xmax[I] << "] " << endl;
  }
  cout << endl;

  cout << "N = " << return_t() << endl;
  cout << setprecision(17)
       << "Vi = " << Vi << ",  Vf = " << return_V()
       << setprecision(6) << endl;
}

void StocDeltaN::sample_plot()
{
  string filename = "sample_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (dim == 1) {
    g.xlabel(string("$N$"));
    g.ylabel(string("$\\phi$"));
    g.plot(Ntraj,x1traj,1,string("b"));
    g.save(filename);
    g.save(filename);
    g.show();
  } else if (dim == 2) {
    g.xlabel(string("$\\phi^1$"));
    g.ylabel(string("$\\phi^2$"));
    g.plot(x1traj,x2traj,1,string("b"));
    g.save(filename);
    g.show();
  }
  g.close();
}

void StocDeltaN::sample_logplot()
{
  string filename = "sample_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (dim == 1) {
    vector<double> x1abs;
    for (auto& x1 : x1traj) {
      x1abs.push_back(fabs(x1));
    }
    
    g.xlabel(string("$N$"));
    g.ylabel(string("$|\\phi|$"));
    g.ylog();
    g.plot(Ntraj,x1abs,1,string("b"));
    g.save(filename);
    g.show();
  } else if (dim == 2) {
    vector<double> x2abs;
    for (auto& x2 : x2traj) {
      x2abs.push_back(fabs(x2));
    }
    
    g.xlabel(string("$\\phi^1$"));
    g.ylabel(string("$|\\phi^2|$"));
    g.ylog();
    g.plot(x1traj,x2abs,1,string("b"));
    g.save(filename);
    g.show();
  }
  g.close();
}

void StocDeltaN::sample_loglinearplot()
{
  string filename = "sample_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (dim == 2) {
    vector<double> x1abs;
    for (auto& x1 : x1traj) {
      x1abs.push_back(fabs(x1));
    }
    
    g.xlabel(string("$|\\phi^1|$"));
    g.ylabel(string("$\\phi^2$"));
    g.xlog();
    g.plot(x1abs,x2traj,1,string("b"));
    g.save(filename);
    g.show();
  }
  g.close();
}

void StocDeltaN::sample_loglogplot()
{
  string filename = "sample_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (dim == 2) {
    g.xlabel(string("$\\phi^1$"));
    g.ylabel(string("$\\phi^2$"));
    g.xlog();
    g.ylog();
    g.plot(x1traj,x2traj,1,string("b"));
    g.save(filename);
    g.show();
  }
  g.close();
}

void StocDeltaN::f1_plot()
{
  string filename = "N_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (dim == 1) {
    g.xlabel(string("$\\phi$"));
    g.ylabel(string("$<N>$"));
    g.plot(site[0],f1,1,string("b"));
    g.save(filename);
    g.show();
  } else if (dim == 2) {
    g.xlabel(string("$\\phi^1$"));
    g.ylabel(string("$\\phi^2$"));
    g.contourf(site[0],site[1],f1,string("$<N>$"));
    g.plot(x1traj,x2traj,3,string("r"));
    g.save(filename);
    g.show();
  }
  g.close();
}

void StocDeltaN::f1_logplot()
{
  string filename = "N_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (dim == 2) {
    vector<double> x2abs;
    for (auto& x2 : x2traj) {
      x2abs.push_back(fabs(x2));
    }
    
    g.xlabel(string("$\\phi^1$"));
    g.ylabel(string("$|\\phi^2|$"));
    g.ylog();
    g.contourf(site[0],site[1],f1,string("$<N>$"));
    g.plot(x1traj,x2abs,3,string("r"));
    g.save(filename);
    g.show();
  }
  g.close();
}

void StocDeltaN::f1_loglinearplot()
{
  string filename = "N_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (dim == 1) {
    g.xlabel(string("$\\phi$"));
    g.ylabel(string("$<N>$"));
    g.xlog();
    g.plot(site[0],f1,1,string("b"));
    g.save(filename);
    g.show();
  } else if (dim == 2) {
    vector<double> x1abs;
    for (auto& x1 : x1traj) {
      x1abs.push_back(fabs(x1));
    }
    
    g.xlabel(string("$|\\phi^1|$"));
    g.ylabel(string("$\\phi^2$"));
    g.xlog();
    g.contourf(site[0],site[1],f1,string("$<N>$"));
    g.plot(x1abs,x2traj,3,string("r"));
    g.save(filename);
    g.show();
  }
  g.close();
}

void StocDeltaN::f1_loglogplot()
{
  string filename = "N_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (dim == 2) {
    vector<double> x1abs, x2abs;
    for (auto& x1 : x1traj) {
      x1abs.push_back(fabs(x1));
    }
    for (auto& x2 : x2traj) {
      x2abs.push_back(fabs(x2));
    }
    
    g.xlabel(string("$|\\phi^1|$"));
    g.ylabel(string("$|\\phi^2|$"));
    g.xlog();
    g.ylog();
    g.contourf(site[0],site[1],f1,string("$<N>$"));
    g.plot(x1abs,x2abs,3,string("r"));
    g.save(filename);
    g.show();
  }
  g.close();
}

void StocDeltaN::g2_plot()
{
  string filename = "dN2_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (dim == 1) {
    g.xlabel(string("$\\phi$"));
    g.ylabel(string("$<\\delta N^2>$"));
    g.ylog();
    g.plot(site[0],g2,1,string("b"));
    g.save(filename);
    g.show();
  } else if (dim == 2) {
    g.xlabel(string("$\\phi^1$"));
    g.ylabel(string("$\\phi^2$"));
    g.log_contourf(site[0],site[1],g2,string("$\\mathrm{log}_{10}<\\delta N^2>$"));
    g.plot(x1traj,x2traj,3,string("r"));
    g.save(filename);
    g.show();
  }
  g.close();
}

void StocDeltaN::g2_logplot()
{
  string filename = "dN2_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (dim == 2) {
    vector<double> x2abs;
    for (auto& x2 : x2traj) {
      x2abs.push_back(fabs(x2));
    }
    
    g.xlabel(string("$\\phi^1$"));
    g.ylabel(string("$|\\phi^2|$"));
    g.ylog();
    g.log_contourf(site[0],site[1],g2,string("$\\mathrm{log}_{10}<\\delta N^2>$"));
    g.plot(x1traj,x2abs,3,string("r"));
    g.save(filename);
    g.show();
  }
  g.close();
}

void StocDeltaN::g2_loglinearplot()
{
  string filename = "dN2_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (dim == 1) {
    g.xlabel(string("$\\phi$"));
    g.ylabel(string("$<\\delta N^2>$"));
    g.xlog();
    g.ylog();
    g.plot(site[0],g2,1,string("b"));
    g.save(filename);
    g.show();
  } else if (dim == 2) {
    vector<double> x1abs;
    for (auto& x1 : x1traj) {
      x1abs.push_back(fabs(x1));
    }
    
    g.xlabel(string("$|\\phi^1|$"));
    g.ylabel(string("$\\phi^2$"));
    g.xlog();
    g.log_contourf(site[0],site[1],g2,string("$\\mathrm{log}_{10}<\\delta N^2>$"));
    g.plot(x1abs,x2traj,3,string("r"));
    g.save(filename);
    g.show();
  }
  g.close();
}

void StocDeltaN::g2_loglogplot()
{
  string filename = "dN2_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (dim == 2) {
    vector<double> x1abs, x2abs;
    for (auto& x1 : x1traj) {
      x1abs.push_back(fabs(x1));
    }
    for (auto& x2 : x2traj) {
      x2abs.push_back(fabs(x2));
    }
    
    g.xlabel(string("$|\\phi^1|$"));
    g.ylabel(string("$|\\phi^2|$"));
    g.xlog();
    g.ylog();
    g.log_contourf(site[0],site[1],g2,string("$\\mathrm{log}_{10}<\\delta N^2>$"));
    g.plot(x1abs,x2abs,3,string("r"));
    g.save(filename);
    g.show();
  }
  g.close();
}

void StocDeltaN::calP_plot()
{
  string filename = "calP_" + model + ".pdf";
  matplotlibcpp g;
  g.open();
  g.xlabel(string("$<N>$"));
  g.ylabel(string("$\\mathcal{P}_\\zeta$"));
  g.ylog();
  g.plot(Ndata,calPdata,1,string("b"));
  g.save(filename);
  g.show();
  g.close();
}

double StocDeltaN::return_intf1()
{
  return Interpolation_f(x,f1);
}

double StocDeltaN::return_intg2()
{
  return Interpolation_f(x,g2);
}
