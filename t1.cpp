#include<iostream>
#include<fstream>
#include <vector>
#include <grace_np.h>
#include "compton.h"
#include "block.h"
using namespace std;

int mcd(int i,int j)
//maximum common divider
{
    if (!i || !j) return(0);
    if (abs(i)==1 || abs(j)==1) return(1);
    int k,i1,j1;
    if (i<j) {
        i1=j;
        j1=i;
    }else{
        i1=i;
        j1=j;
    }
    do{
        k=i1%j1;
        i1=j1;
        j1=k;
    }while (k);
    return(i1);
}

int main()
{
    vector<unsigned int> dist_theta;
    int theta_steps=1000;
    dist_theta.resize(theta_steps+1);
    double dtheta=theta_steps/M_PI;
    for (int i=0;i<=theta_steps;i++) dist_theta.at(i)=0;
    double t0;
    int j;
    if (GraceOpenVA("xmgrace", 16384, "-nosafe", "-noask","-nosigcatch", NULL)==-1){
        cerr<<"Can't run Grace. \n";
        exit(1);
    }
    GracePrintf("s0 on");
    GracePrintf("s0 type xy");
    GracePrintf("s1 on");
    GracePrintf("s1 type xy");
    GracePrintf("s1 hidden true");
    GracePrintf("s0 symbol 1");
    GracePrintf("s0 symbol size 0.3");
    GracePrintf("s0 line linestyle 0");
    GracePrintf("world xmin 0");
    GracePrintf("world xmax %g",(double) M_PI);
    GracePrintf("xaxis  tick spec type both");
    GracePrintf("xaxis  tick spec 10");
    unsigned int xaxisDiv=8;
    for (unsigned int ii=0;ii<=xaxisDiv;ii++){
        if ( ii & (unsigned int)1) {
            GracePrintf("xaxis tick minor %d,%g",ii,(double) ii/xaxisDiv*M_PI);
        }else {
            GracePrintf("xaxis tick major %d,%g",ii,(double) ii/xaxisDiv*M_PI);
            if (!ii){
                GracePrintf("xaxis ticklabel %d,\"0\"",ii);
                continue;
            }
            int ii1=mcd(ii,xaxisDiv);
            //cout<<ii<<' '<<xaxisDiv<<' '<<ii1<<endl;
            if (xaxisDiv==ii){
                GracePrintf("xaxis ticklabel %d,\"\\f{Symbol}p\"",ii);
                continue;
            }
            GracePrintf("xaxis ticklabel %d,\"%d/%d\\f{Symbol}p\"",ii,ii/ii1,xaxisDiv/ii1);
        }
    }

    int l1=1000,l2=50000;
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
                if (fabs(p0.o.get_phi())<5e-3) dist_theta.at( (int)(p0.o.get_theta() * dtheta+ 0.5))++;
            }
        }
        cout<<"ii="<<ii<<endl;
        GracePrintf("kill s0");
        GracePrintf("s0 on");
        GracePrintf("s0 type xy");
        GracePrintf("s0 symbol 1");
        GracePrintf("s0 symbol size 0.3");
        GracePrintf("s0 line linestyle 0");
        double fac2=fac1/ii;
        ofstream out1("t1.out");
        for (int i=1;i<=theta_steps;i++){
            double xi=i*M_PI/theta_steps,yi=(double) dist_theta.at(i)*fac2;
            out1<<xi<<' '<<yi<<endl;
            if (dist_theta.at(i)) GracePrintf("s0 point %g,%g",xi,yi);
        }
        out1.close();
        GracePrintf("autoscale yaxes");
        GracePrintf("autoticks");
        GracePrintf("redraw");
    }





    return(0);
}
