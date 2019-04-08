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
  vector< vector< vector<double> > > site, hI;
  // site[xp][I][index] : phase-space value,  hI[xp][I][index] : step size
  vector< vector<int> > siteNo; // siteNo[xp][I] : site number in xpI
  vector<double> *ff, *f_next; // solved function
  
  // ---------------------------------------------------------
  /*
  num_p[num][xp][I] : site number positively next to num in xpI direction
  num_pp[num][xpI][I][xpJ][J] : positively next to positively next in xpI and xpJ direction
  hp[num][xp][I] : site step positively next to num in xpI direction
  set in BoundaryCondition
  */
  vector< vector< vector<int> > > num_p, num_m;
  vector< vector< vector< vector< vector<int> > > > > num_pp, num_pm, num_mm;
  vector< vector< vector<double> > > hp, hm;
  // ---------------------------------------------------------
  
  double maxstep,tol,rhoc,xpdim,Idim;
  int volume, funcNo;
  vector<int> xpvol;
  vector<bool> Omega;

public:
  JacobiPDE(){}
  // Params[0] = maxstep, Params[1] = tol, Params[2] = funcNo, Params[3] = rhoc
  JacobiPDE(vector< vector< vector<double> > > &Site, vector<double> &Params);
  double PDE_1step(int num, int func);
  void PDE_solve(int func);
  int Ind2No(vector< vector<int> > &index); // index set to site number
  int No2Ind(int num, int xp, int I); // site number to index of xpI
  double No2PSV(int num, int xp, int I); // site number to phase-space value of xpI
  int ceilXP(int xp, int I, vector< vector<double> > &psv);
  double Interpolation_f(vector< vector<double> > &psv, int func);
  void export_fg(string filename);
  virtual double H(vector<double> &X, vector<double> &P);
  virtual double V(vector<double> &X);
  virtual double VI(vector<double> &X, int I);
  virtual double metric(vector<double> &X, int I, int J);
  virtual double inversemetric(vector<double> &X, int I, int J);
  virtual double affine(vector<double> &X, int I, int J, int K);
  virtual double derGamma(vector<double> &X, int I, int J, int K, int L); // Gamma^I_{JK,L}
  virtual double DI(int xp, int I, vector< vector<double> > &psv);
  virtual double DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv);
  virtual double CC(int num, vector< vector<double> > &psv, int func);
  virtual void BoundaryCondition();
  virtual bool EndSurface(vector< vector<double> > &psv);
};

#endif
