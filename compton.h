#pragma once
#include<iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#define int_steps 4000
#define E_electron 510.998918
using namespace std;
double e_new(double,void*),dSigma_dOmega(double, void*);
int mcd(int,int);

class thetaDistribution //photon energy
{
private:
    vector<double> thetai,fi;
    gsl_interp_accel *acc;
    gsl_spline *spline;

public:
    thetaDistribution(double);
    double theta(void);
    double e0,sigma_total;
};
