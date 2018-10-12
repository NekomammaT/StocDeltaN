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
  double t, t0,
    A0[3][3],A1[3][3],B0[3][3],B1[3][3],Alpha[3],Beta1[3],Beta2[3],C0[3],C1[3];
  vector<double> x,xi,p,pi,dW,ax[3],ap[3],H0x[3],H0p[3];
  vector< vector<double> > bx[3],bp[3],Hkx[3],Hkp[3],aIs,vIs;

public:
  SRKintegrater(){}
  SRKintegrater(vector<double> &Xi, vector<double> &Pi, double T0, int NoiseDim);
  void SRK2(double dt);
  void coeff(double dt, int step);
  double e1(vector<double> &X, vector<double> &P);
  double return_H();
  double return_t();
  double return_phi(int I);
  double return_pi(int I);
  double return_e1();
  double vielbein(vector<double> &X, vector<double> &P, int I, int alpha);
  double eIsigma(vector<double> &X, vector<double> &P, int I);
  double eIs(vector<double> &X, vector<double> &P, int I, int alpha);
  virtual double H(vector<double> &X, vector<double> &P);
  virtual double V(vector<double> &X);
  virtual double VI(vector<double> &X, int I);
  virtual double metric(vector<double> &X, int I, int J);
  virtual double inversemetric(vector<double> &X, int I, int J);
  virtual double affine(vector<double> &X, int I, int J, int K);
  virtual double Dphi(vector<double> &X, vector<double> &P, int I);
  virtual double Dpi(vector<double> &X, vector<double> &P, int I);
  virtual double PhiNoise(vector<double> &X, vector<double> &P, int I, int alpha);
  virtual double PiNoise(vector<double> &X, vector<double> &P, int I, int alpha);
};

double Uniform();
double rand_normal(double mu, double sigma);

#endif
