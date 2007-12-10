//scattering in reflectivity geometry
//background geometry, beam directly passing the cylinder sample
#include<iostream>
#include<fstream>
#include <vector>
#include <grace_np.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "compton.h"
#include "block.h"
#include "randomInit.h" //init random number generator
using namespace std;

int mcd(int i,int j)
//maximum common divider
{
    char randomBuffer[257];
    randomInit(randomBuffer); //init random number generator
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
        double s3h=6./640.,s3v=6./640.;
    vector<unsigned int> dist_theta;
    int theta_steps=(int)(M_PI/s3h+0.5);
    dist_theta.resize(theta_steps+2);
    for (unsigned int i=0;i<dist_theta.size();i++) dist_theta.at(i)=0;
    double dtheta=theta_steps/M_PI;
    int j;
    /*
    unsigned int xaxisDiv=8;
    for (unsigned int ii=0;ii<=xaxisDiv;ii++){
        if ( ii & (unsigned int)1) {
        }else {
            if (!ii){
                continue;
            }
            int ii1=mcd(ii,xaxisDiv);
            //cout<<ii<<' '<<xaxisDiv<<' '<<ii1<<endl;
            if (xaxisDiv==ii){
                continue;
            }
        }
    }

    */
    string fn("si-background0.txt");
    cout<<"Deleting output file "<<fn<<endl;
    unlink(fn.c_str());
    int l1=100,l2=50000000;
    ofstream out1;
    double l3=0.;
    int ii=1;
    struct rusage r_start,r_end;
                static thetaDistribution tP0(30.);
            static photon p0;
            static double dthetaPriv=dtheta;

    getrusage(RUSAGE_SELF, &r_start);
    while (ii<l1){
        //double t1=ii*fac1;
//#pragma omp parallel for default(shared) private( j,t0) schedule(static)
        for (int i=0;i<l2;i++){
            if (p0.init()) continue;
            while (( j=p0.propagate(tP0.theta()))==1);
            // if( i && ((i>>18 ) <<18) == i) cout<<"i="<<i<<endl;
            if (j==-1 ) {//scattered out of sample
                //    cout<<i<<' '<<p0.o.theta<<endl;
                if (fabs(p0.o.get_phi())<s3v) dist_theta.at( (int)(p0.o.get_theta() * dthetaPriv+ 0.5))++;
            }
        }
        getrusage(RUSAGE_SELF, &r_end); //get running time
        l3=(double) ii*l2;
        cout<<ii<<' '<<dist_theta.at(theta_steps>>2)<<' '<<l3<<' '<<l3/(r_end.ru_utime.tv_sec -r_start.ru_utime.tv_sec+ double(1e-6)*(r_end.ru_utime.tv_usec -r_start.ru_utime.tv_usec ) )<<" P/s\n";
        //r_start=r_end;
        if (dist_theta.at(theta_steps>>2)>6000) {
            out1.open(fn.c_str(),fstream::app);
            dtheta=1./dtheta;
            double fac2=dtheta/(s3v*s3h*l3);
            for (int i1=1;i1<=theta_steps;i1++){
                double xi=i1*dtheta,yi=(double) dist_theta.at(i1)*fac2;
                out1<<xi<<' '<<yi<<' '<< dist_theta.at(i1)<<' '<<l3<<endl;
            }
            out1.close();
            break;
        }
        ii++;
    }





    return(0);
}
