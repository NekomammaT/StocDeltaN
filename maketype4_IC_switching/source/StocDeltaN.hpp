// stochastic-delta N solver

#ifndef INCLUDED_StocDeltaN_hpp_
#define INCLUDED_StocDeltaN_hpp_

#include "JacobiPDE.hpp"
#include "SRK32.hpp"

class StocDeltaN: virtual public JacobiPDE, virtual public SRKintegrater
{
protected:
  string model; // model name
  double timestep,Nmax,deltaN;
  // timestep : time step dN for SRK
  // Nmax : max <N> for which power spectrum is calculated
  // deltaN : <N>'s step for power spectrum calculation
  int recursion,xpdim,Idim;
  // recursion : number of sample paths for power spectrum calculation
  // xpdim : xpdim=1 -> slow-roll field-space, xpdim=2 -> full phase-space
  // Idim : number of inflatons
  vector<double> x1traj,p1traj,x2traj,Ntraj,Ndata,calPdata;
  // x1traj : trajectory of phi^1 for plotting
  // p1traj : trajectory of pi_1
  // x2traj : trajectory of phi^2
  // Ntraj : trajectory of e-folding time N
  // Ndata : data variable for <N> in power spectrum calculation
  // calPdata : data variable for power spectrum
  
public:
  StocDeltaN(){}
  StocDeltaN(string Model, vector< vector< vector<double> > > &Site,
	     vector< vector<double> > &XPi, double T0, vector<double> &Params);
  // Params[0] = maxstep, Params[1] = tol, Params[2] = funcNo, Params[3] = rhoc,
  // Params[4] = NoiseDim, Params[5] = timestep, Params[6] = Nmax, Params[7] = deltaN,
  // Params[8] = recursion
  void solve(); // execute stochastic-delta N
  void sample(); // obtain 1 sample path

  // -------------- plot sample path -------------
  // Idim=1 & xpdim=1 : N vs phi
  // Idim=1 & xpdim=2 : phi vs pi
  // Idim=2           : phi vs psi
  void sample_plot();
  void sample_logplot();
  void sample_loglinearplot();
  void sample_loglogplot();
  // ---------------------------------------------

  // ------------- plot <N> or <delta N^2> -------------
  // xpdim=1 & Idim=1 : phi vs <N> or <delta N^2>
  // xpdim=1 & Idim=2 : contour of <N> or log_10 <delta N^2> in (phi,psi) plane
  // xpdim=2 & Idim=1 : contour of <N> or log_10 <delta N^2> in (phi,pi) plane
  void f_plot(int func);
  void f_logplot(int func);
  void f_loglinearplot(int func);
  void f_loglogplot(int func);
  // ---------------------------------------------------
  
  void calP_plot(); // plot power spectrum of curvature perturbation
  double return_intf(int func); // return interpolated <N> or <delta N^2> at the current point of sample path

  virtual double V(vector<double> &X); // potential
  virtual double VI(vector<double> &X, int I); // \partial_I V
  virtual double metric(vector<double> &X, int I, int J); // field-space metric G_IJ
  virtual double inversemetric(vector<double> &X, int I, int J); // inverse field-space metric G^IJ
  virtual double affine(vector<double> &X, int I, int J, int K); // Christoffel symbol Gamma^I_JK
  virtual double derGamma(vector<double> &X, int I, int J, int K, int L); // Gamma^I_{JK,L}

  // ----------------- validate to modify them if want ---------------
  //virtual void BoundaryCondition(double Ncut);
  //virtual double H(vector<double> &X, vector<double> &P);
  //virtual double DI(int xp, int I, vector< vector<double> > &psv);
  //virtual double DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv);
  //virtual double gIa(int xp, int I, int alpha, vector< vector<double> > &psv);
  //virtual double CC(int num, vector< vector<double> > &psv, int func);
  //virtual bool EndSurface(vector< vector<double> > &psv);
  // -----------------------------------------------------------------
};

#endif
