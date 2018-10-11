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
