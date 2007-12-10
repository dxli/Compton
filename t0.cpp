#include<iostream>
#include "compton.h"
#include "block.h"
using namespace std;

int main()
{
    thetaDistribution tP0(30.);
//cout<<tP0.sigma_total<<endl;
    vector<unsigned int> dist;
    int llength=200;
    double dx=M_PI/llength;
    dist.resize(llength+1);
    int steps=2000000;
    for (int i=0;i<steps;i++){
        dist.at((int)(tP0.theta()/dx))++;
    }
    for (int i=0;i<llength;i++){
        cout<<i*dx<<' '<<(double) dist.at(i)/steps*llength/M_PI<<' '<<dSigma_dOmega(i*dx,&(tP0.e0))*2*M_PI/tP0.sigma_total<<endl;
    }
    return(0);
}
