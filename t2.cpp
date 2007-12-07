#include<iostream>
#include<fstream>
#include <vector>
#include <grace_np.h>
#include "compton.h"
#include "block.h"
using namespace std;

int main()
{
    vector<unsigned int> dist_theta;
    int theta_steps=1000;
    dist_theta.resize(theta_steps+1);
    double dtheta=theta_steps/M_PI;
    for (int i=0;i<=theta_steps;i++) dist_theta.at(i)=0;
    double t0;
    int j;
    int l1=1000,l2=10000000;
    double fac1=1./l2;
    for (int ii=1;ii<l1;ii++){
#pragma omp parallel for private( j,t0 ) schedule(static)
        for (int i=0;i<l2;i++){
            static thetaDistribution tP0(30.);
            static photon p0;
            if (p0.init()) continue;
            do{
                t0=tP0.theta();
                j=p0.propagate(t0);
                // cout<<p0.r.y<<endl;
            }while (j==1);
            // if( i && ((i>>18 ) <<18) == i) cout<<"i="<<i<<endl;
            if (j==-1 ) {//scattered out of sample
                //    cout<<i<<' '<<p0.o.theta<<endl;
                if (fabs(p0.o.phi)<5e-3) dist_theta.at( (int)(p0.o.theta * dtheta+ 0.5))++;
            }
        }
        cout<<"ii="<<ii<<endl;
        ofstream out1("t2.out");
        double fac2=fac1/ii;
        for (int i=1;i<=theta_steps;i++){
            if (dist_theta.at(i)) out1<<i*M_PI/theta_steps<<' '<<(double) dist_theta.at(i)*fac2<<endl;
        }
        out1.close();
    }





    return(0);
}
