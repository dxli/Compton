//scattering in pass-through geometry
//pass through geometry, beam pass through Si crystal
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
    vector<double> dist_theta;
    int theta_steps=(int)(M_PI/s3h+0.5);
    double itheta_steps=theta_steps/M_PI;

    dist_theta.resize(theta_steps+2);
    for (unsigned int i=0;i<dist_theta.size();i++) dist_theta.at(i)=0.;
    int j;
    string fn("si-Compton-pass-new.txt");
    cout<<"Deleting output file "<<fn<<endl;
    unlink(fn.c_str());
    int l1=100,l2=500000000;
    ofstream out1;
    double l3=0.;
    int ii=1;
    double      muC= 0.15 * 2.33; // in cm^-1, compton
    double  muTotal=1.31 * 2.33; // in cm^-1, total
    double compton_ratio_factor=muC/muTotal; //weight factor for each scattering
    double psum;
    double pEn=30.; //KeV
    //     double pLambda=E_to_l(pEn);
    //   double pk0=2*M_PI/pLambda;
    //double fac1=1./(2.*pk0*l1); // q_z range up to 1 angstrom^-1
    thetaDistribution tP0(pEn);
    photon p0;

    struct rusage r_start,r_end;
    getrusage(RUSAGE_SELF, &r_start);
    psum=0;
    while (ii<l1){
        //double t1= -ii*fac1;
//#pragma omp parallel for private( j,t0 ) reduction(+:sum) schedule(static)
        for (int i=0;i<l2;i++){
            //if (p0.initRefBuried(t1)) continue;
            p0.initPass();

            while ((j=p0.propagatePass(tP0.theta()))==1);
            // if( i && ((i>>18 ) <<18) == i) cout<<"i="<<i<<endl;
            if (j==-1 ) {//scattered out of sample
                //    cout<<i<<' '<<p0.o.theta<<endl;
                // cout<<i<<' '<<tP0.theta()<<' '<<p0.o.st<<' '<<p0.o.sp<<endl;
                if (p0.scattered && fabs(p0.o.sp)<s3v) dist_theta.at((unsigned int) (p0.o.get_theta()*itheta_steps+0.5)) += pow(compton_ratio_factor,(int) p0.scattered);
            }
        }
        getrusage(RUSAGE_SELF, &r_end); //get running time
        l3+=l2;
        cout<<ii<<' '<<dist_theta.at(ii)<<' '<<l3<<' '<<l2/(r_end.ru_utime.tv_sec -r_start.ru_utime.tv_sec+ double(1e-6)*(r_end.ru_utime.tv_usec -r_start.ru_utime.tv_usec ) )<<" P/s\n";
        r_start=r_end;
        /*
        */
        out1.open(fn.c_str());
        double fac2=1./(s3v*s3h*l3);
        double dx=1./itheta_steps;
        double x0=0.5*dx;
        for (unsigned int i2=0;x0<M_PI;i2++){
            out1<<x0<<' '<<dist_theta.at(i2)*fac2<<' '<<dist_theta.at(ii)<<' '<<l3<<endl;
            x0+=dx;
        }
        out1<<dist_theta.size()/itheta_steps<<' '<<dist_theta.at(ii)*fac2<<' '<<s3v<<' '<<s3h<<endl;
        out1.close();
        ii++;
    }
    return(0);
}
