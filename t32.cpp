//scattering in reflectivity geometry
//reflectivity geometry, beam pass through Si crystal
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

int main()
{
            double s3h=20./640.,s3v=2./640.;// slit size
            vector<unsigned int> dist_theta;
            int theta_steps=(int)(M_PI/s3h+0.5);

    dist_theta.resize(theta_steps+2);
    for (unsigned int i=0;i<dist_theta.size();i++) dist_theta.at(i)=0;
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
    string fn("si-Compton-RefBuried-new.txt");
    cout<<"Deleting output file "<<fn<<endl;
    unlink(fn.c_str());
    int l1=100,l2=500000000;
    ofstream out1;
    double l3=0.;
    int ii=1;
        unsigned int sum;
        double pEn=30.; //KeV
        double pLambda=E_to_l(pEn);
        double pk0=2*M_PI/pLambda;
    double fac1=1./(2.*pk0*l1); // q_z range up to 1 angstrom^-1
            thetaDistribution tP0(pEn);
            photon p0;

    struct rusage r_start,r_end;
    getrusage(RUSAGE_SELF, &r_start);
    while (ii<l1){
        double t1= -ii*fac1;
        sum=0;
//#pragma omp parallel for private( j,t0 ) reduction(+:sum) schedule(static)
        for (int i=0;i<l2;i++){
            //if (p0.initRefBuried(t1)) continue;
            p0.initRefBuried(t1);

                while((j=p0.propagateRef(tP0.theta()))==1);
            // if( i && ((i>>18 ) <<18) == i) cout<<"i="<<i<<endl;
            if (j==-1 ) {//scattered out of sample
                //    cout<<i<<' '<<p0.o.theta<<endl;
                   // cout<<i<<' '<<tP0.theta()<<' '<<p0.o.st<<' '<<p0.o.sp<<endl;
                if (p0.scattered && fabs(p0.o.sp+t1)<s3v && fabs(p0.o.st)<s3h) sum++;
            }
        }
        dist_theta.at(ii)+=sum;
        getrusage(RUSAGE_SELF, &r_end); //get running time
        l3+=l2;
        cout<<ii<<' '<<dist_theta.at(ii)<<' '<<l3<<' '<<l2/(r_end.ru_utime.tv_sec -r_start.ru_utime.tv_sec+ double(1e-6)*(r_end.ru_utime.tv_usec -r_start.ru_utime.tv_usec ) )<<" P/s\n";
        r_start=r_end;
        if (dist_theta.at(ii)>1000) {
            /*
            */
            out1.open(fn.c_str(),fstream::app);
            double fac2=1./(cos(t1)*s3v*s3h*l3);
            double xi=30/12.4*4*M_PI*sin(t1),yi=(double) dist_theta.at(ii)*fac2;
            out1<<ii<<' '<<dist_theta.at(ii)<<' '<<l3<<' '<<xi<<' '<<yi<<endl;
            ii++;
            out1.close();
            l3=0.;
        }
    }





    return(0);
}
