#ifndef INCLUDED_JacobiPDE_hpp_
#define INCLUDED_JacobiPDE_hpp_

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

class JacobiPDE
{
protected:
  vector<double> f1, f1_next, g2, g2_next; // argument: number
  vector< vector<double> > site[2], hI[2]; // argument: x or p, I, siteNo
  vector<double> FPoint[2]; // argument: x or p, I
  vector<bool> Omega; // argument: number
  vector<int> ind[2]; // argument: x or p, I
  vector< vector<int> > numXp, numXm, numPp, numPm; // argument: number, I
  vector< vector< vector<int> > > numXXpp, numXXpm, numXXmm,
    numXPpp, numXPpm, numXPmp, numXPmm, numPPpp, numPPpm, numPPmm; // argument: number, I, J
  vector< vector<double> > hp[2], hm[2]; // argument: x or p, number, I
  vector<int> siteNo[2]; // argument: x or p, I
  double rhoc;
  int dim, volume, pvol;
  
public:
  JacobiPDE(){}
  JacobiPDE(vector< vector<double> > Site[], double Rhoc);
  double PDE(int number, int n);
  void PDE_solve(int maxstep, double tol, int n);
  int Ind2No(vector<int> index[]);
  int No2XInd(int number, int I);
  int No2PInd(int number, int I);
  double No2X(int number, int I);
  double No2P(int number, int I);
  int ceilX(vector<double> &X, int I);
  int ceilP(vector<double> &X, int I);
  double Interpolation_f(vector<double> &X, vector<double> &P, vector<double> &f);
  double return_f1(vector<int> index[]);
  double return_g2(vector<int> index[]);
  void export_fg(string filename);
  virtual double H(vector<double> &X, vector<double> &P);
  virtual double V(vector<double> &X);
  virtual double VI(vector<double> &X, int I);
  virtual double metric(vector<double> &X, int I, int J);
  virtual double inversemetric(vector<double> &X, int I, int J);
  virtual double affine(vector<double> &X, int I, int J, int K);
  virtual double Dphi(vector<double> &X, vector<double> &P, int I);
  virtual double Dpi(vector<double> &X, vector<double> &P, int I);
  virtual double Dphiphi(vector<double> &X, vector<double> &P, int I, int J);
  virtual double Dphipi(vector<double> &X, vector<double> &P, int I, int J);
  virtual double Dpipi(vector<double> &X, vector<double> &P, int I, int J);
};

#endif
