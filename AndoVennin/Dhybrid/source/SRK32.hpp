// 3-step strong order 2 stochastic Runge-Kutta solver

#ifndef INCLUDED_SRK32_hpp_
#define INCLUDED_SRK32_hpp_

#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <time.h>

using namespace std;

class SRKintegrater
{
protected:
  double t,t0, // time t and its initial value t0
    A0[3][3],A1[3][3],B0[3][3],B1[3][3],Alpha[3],Beta1[3],Beta2[3],C0[3],C1[3]; // params. for stochatic Runge-Kutta
  vector<double> dW; // dW[alpha], alpha : noise d.o.f.
  vector< vector<double> > xx,xxi,aa[3],H0[3]; // xx[xp][I], xp: phi or pi, I : field d.o.f.
  // xx, xxi : phase-space value and its initial condition
  // aa, H0 : for SRK
  vector< vector<double> > aIs,vIs; // aIs[I][alpha] only for Idim == noisedim
  // for calc. of adiabatic direction
  vector< vector< vector<double> > > bb[3],Hk[3]; // bb[3][alpha][xp][I]
  // for SRK 
  int xpdim,Idim,noisedim; // xpdim=1 : slow roll field-space, xpdim=2 : full phase-space
  // Idim : number of inflatons, noisedim : noise d.o.f.

public:
  SRKintegrater(){}
  SRKintegrater(vector< vector<double> > &XPi, double T0, int NoiseDim);
  // XPi : initial phase space value, T0 : initial time, NoiseDim : noise d.o.f.
  void SRK2(double dt); // execute 1 SRK step with time step dt
  void coeff(double dt, int step); // set coefficients of SDE
  void init_txp(); // initialize t and phase-space value to initial values
  void set_txp(double T, vector< vector<double> > &PSV); // set t and phase-space value by arbitrary values T and PSV
  double e1(vector<double> &X, vector<double> &P); // first SR param. -\dot{H}/H^2
  double return_H(); // return current Hubble
  double return_V(); // return current potential
  double return_t(); // return current time
  double return_xp(int xp, int I); // return current phase-space value in xpI direction
  double return_e1(); // return current -\dot{H}/H^2
  double vielbein(vector< vector<double> > &XP, int I, int alpha); // projection from field coordinate xpI to noise frame alpha
  double eIsigma(vector< vector<double> > &XP, int I); // projection to adiabatic direction
  double eIs(vector< vector<double> > &XP, int I, int alpha); // projection to entropic direction labeled by alpha
  double eta_perp(vector< vector<double> > &XP); // turning param. \eta_\perp = -V_s/(H\dot{\phi}). Only for 2-field model
  double return_etaperp(); // return current \eta_\perp
  virtual double H(vector<double> &X, vector<double> &P); // Hubble parameter
  virtual double V(vector<double> &X); // potential
  virtual double VI(vector<double> &X, int I); // \partial_I V
  virtual double metric(vector<double> &X, int I, int J); // field-space metric G_IJ
  virtual double inversemetric(vector<double> &X, int I, int J); // inverse field-space metric G^IJ
  virtual double affine(vector<double> &X, int I, int J, int K); // Christoffel symbol Gamma^I_JK
  virtual double derGamma(vector<double> &X, int I, int J, int K, int L); // Gamma^I_{JK,L}

  // ---------------------------------------------------------------
  virtual double DI(int xp, int I, vector< vector<double> > &psv);
  virtual double DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv);
  virtual double gIa(int xp, int I, int alpha, vector< vector<double> > &psv);
  // solve SDE : dXI = DI dN + gIa dWa
  // diffusion g^I_alpha g^J_alpha = D^IJ
  // ---------------------------------------------------------------
};

double Uniform(); // uniform random number (0,1)
double rand_normal(double mu, double sigma);
// normal distribution
// \exp(-(x-\mu)^2/2/\sigma^2)/\sqrt{2\pi\sigma^2}

#endif
