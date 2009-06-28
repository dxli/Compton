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
#define g_l3 550. //g_l3 in mm

int main()
{
    double s3h=20./(g_l3),s3v=2./(g_l3), // slit size
                            qz_max=2.5; // in angstrom^-1
    int qz_steps=20;
    double dqz=qz_max/qz_steps;
    vector<double> dist_theta;
    char randomBuffer[257];
    randomInit(randomBuffer); //init random number generator

    dist_theta.resize(qz_steps+2);
    for (unsigned int i=0;i<dist_theta.size();i++) dist_theta.at(i)=0.;
    int j;
    string fn("si-Compton-RefBuried-1.txt");
    cout<<"Deleting output file "<<fn<<endl;
    unlink(fn.c_str());
    int l2=5000000000; // photons per loop
    ofstream out1;
    double l3=0.;
    int ii=1;
    double      muC= 0.15 * 2.33; // in cm^-1, compton
    double  muTotal=1.31 * 2.33; // in cm^-1, total
    double compton_ratio_factor=muC/muTotal; //weight factor for each scattering
    double psum;
    double pEn=30.; //KeV
    double pLambda=E_to_l(pEn); //wavelength
    double pk0=2*M_PI/pLambda;
    double fac1=-dqz/(2.*pk0); // dsin\theta
    thetaDistribution tP0(pEn);
    photon p0(0.25); //sample diameter 0.25"
    struct rusage r_start,r_end;
            double fac2=1./(s3v*s3h*l2);
    vector<double> pEp,pEp2;
    pEp.resize(qz_steps+2);
    pEp2.resize(qz_steps+2);
    for(int i=0;i<qz_steps+2;i++){
            pEp.at(i)=0;
            pEp2.at(i)=0;
    }
    for(int i=1;i<=1000000;i++){
    getrusage(RUSAGE_SELF, &r_start);
    psum=0;
    for(ii=1;ii<=qz_steps;ii++){
        double t1= ii*fac1;
//#pragma omp parallel for private( j,t0 ) reduction(+:sum) schedule(static)
        for (int i=0;i<l2;i++){
            p0.initRefBuried(t1);
            while ((j=p0.propagateRef(tP0.theta()))==1);
            if (j==-1 ) {//scattered out of sample
                if (p0.scattered && fabs(fabs(p0.o.sp)+t1)<s3v && fabs(p0.o.st)<s3h) psum += pow(compton_ratio_factor,(int) p0.scattered);
            }
        }
        getrusage(RUSAGE_SELF, &r_end); //get running time
        cout<<ii<<' '<<psum<<' '<<l2<<' '<<l2/(r_end.ru_utime.tv_sec -r_start.ru_utime.tv_sec+ double(1e-6)*(r_end.ru_utime.tv_usec -r_start.ru_utime.tv_usec ) )<<" P/s\n";
        r_start=r_end;
        /*
        */
        //counts for cross-section
        psum *=fac2;
        pEp.at(ii) += psum;
        pEp2.at(ii) += psum*psum;
    }
    out1.open(fn.c_str());
    for(ii=1;ii<=qz_steps;ii++){
            double a=pEp.at(ii)/i;
            out1<<ii*dqz<<' '<<a<<' '<<pEp2.at(ii)/i-a*a<<endl;
    }
    out1.close();
    cout<<"(qz cross-section) Written results to "<<fn<<endl;
    }
    return(0);
}
