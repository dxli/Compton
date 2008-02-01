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
    double s3h=20./640.,s3v=2./640., // slit size
                            qz_max=2.; // in angstrom^-1
    int qz_steps=100;
    double dqz=qz_max/qz_steps;
    vector<double> dist_theta;
    char randomBuffer[257];
    randomInit(randomBuffer); //init random number generator

    dist_theta.resize(qz_steps+2);
    for (unsigned int i=0;i<dist_theta.size();i++) dist_theta.at(i)=0.;
    int j;
    string fn("si-Compton-Ref-bg.txt");
    cout<<"Deleting output file "<<fn<<endl;
    unlink(fn.c_str());
    int l2=500000000; // photons per loop
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
    photon p0(0.125);
    struct rusage r_start,r_end;
    getrusage(RUSAGE_SELF, &r_start);
    psum=0;
        double t1= sin(0.35*M_PI/180.);
    while (ii<qz_steps){
//#pragma omp parallel for private( j,t0 ) reduction(+:sum) schedule(static)
        for (int i=0;i<l2;i++){
            p0.initRef(t1);
            while ((j=p0.propagateRef(tP0.theta()))==1);
            if (j==-1 ) {//scattered out of sample
                if (p0.scattered &&  fabs(p0.o.st)<s3h && p0.o.sp<0. ) {
                        dist_theta.at((int) (fabs(p0.o.sp)*qz_steps))++;
                }
            }
        }
        getrusage(RUSAGE_SELF, &r_end); //get running time
        l3+=l2;
        cout<<ii<<' '<<psum<<' '<<l3<<' '<<l2/(r_end.ru_utime.tv_sec -r_start.ru_utime.tv_sec+ double(1e-6)*(r_end.ru_utime.tv_usec -r_start.ru_utime.tv_usec ) )<<" P/s\n";
        r_start=r_end;
        /*
        */
            out1.open(fn.c_str(),fstream::app);
            double fac2=1./qz_steps,fac3=1./(l3*s3h*fac2);
            for(int i1=0;i1<=qz_steps;i1++){
            out1<<i1*fac2<<' '<<dist_theta.at(i1)*fac3<<endl;
            }
            out1.close();
            l3=0;
    }
    return(0);
}
