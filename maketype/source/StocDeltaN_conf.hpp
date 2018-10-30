#ifndef INCLUDED_StocDeltaN_hpp_
#define INCLUDED_StocDeltaN_hpp_

#include "JacobiPDE_conf.hpp"
#include "SRK32_conf.hpp"

class StocDeltaN: virtual public JacobiPDE, virtual public SRKintegrater
{
protected:
  string model;
  double maxstep, tol, timestep, Nmax, deltaN;
  int recursion;
  vector<double> Ntraj, x1traj, x2traj, Ndata, calPdata;
  
public:
  StocDeltaN(){}
  StocDeltaN(string Model,
	     vector< vector<double> > &Site, double Rhoc,
	     vector<double> &Xi, double T0,
	     double Maxstep, double Tol, int Recursion,
	     double Timestep, double NNmax, double DeltaN);
  void init_fn();
  void init_txp();
  void solve();
  void sample();
  void sample_plot();
  void sample_logplot();
  void sample_loglinearplot();
  void sample_loglogplot();
  void f1_plot();
  void f1_logplot();
  void f1_loglinearplot();
  void f1_loglogplot();
  void g2_plot();
  void g2_logplot();
  void g2_loglinearplot();
  void g2_loglogplot();
  void calP_plot();
  double return_intf1();
  double return_intg2();
  virtual double V(vector<double> &X);
  virtual double VI(vector<double> &X, int I);
  virtual double metric(vector<double> &X, int I, int J);
  virtual double inversemetric(vector<double> &X, int I, int J);
  virtual double affine(vector<double> &X, int I, int J, int K);
  virtual double PhiNoise(vector<double> &X, int I, int alpha);
  virtual double Dphi(vector<double> &X, int I);
  virtual double Dphiphi(vector<double> &X, int I, int J);
};

#endif
