//scattering in pass-through geometry
//pass through geometry, beam pass through Si crystal
#include<iostream>
#include<fstream>
#include<sstream>
#include <vector>
#include <grace_np.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "compton.h"
#include "block.h"
#include "randomInit.h" //init random number generator
using namespace std;

int main(int argc, char *argv[])
{
    char randomBuffer[257];
    randomInit(randomBuffer); //init random number generator
    double pixelHeight=0.172/900.// angular now
            ,s5h=2/900.;// slit size
    vector<double> dist_theta;
    int theta_steps=(int)(1/pixelHeight+0.5);
    double itheta_steps=theta_steps/(1.);
    dist_theta.resize(theta_steps+2);
    for (unsigned int i=0;i<dist_theta.size();i++) dist_theta.at(i)=0.;
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
    double diameter=0.;
    if (argc>=2) { // read sample block diameter from command line
        istringstream is0(argv[1]);
        is0>>diameter;
        if ( !( diameter>=0.01 && diameter < 1000.)) diameter=0.25;
    } else diameter=0.25; // all units in inches
    ostringstream os0;
    os0<<"si-gixos-"<<diameter<<".txt";
    string fn(os0.str().c_str());
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
    photon p0(diameter);

    struct rusage r_start,r_end;
    getrusage(RUSAGE_SELF, &r_start);
    psum=0;
    while (ii<l1){
        //double t1= -ii*fac1;
//#pragma omp parallel for private( j,t0 ) reduction(+:sum) schedule(static)
        for (int i=0;i<l2;i++){
            //if (p0.initRefBuried(t1)) continue;
            p0.initGixos();

            while ((j=p0.propagateRef(tP0.theta()))==1);
            // if( i && ((i>>18 ) <<18) == i) cout<<"i="<<i<<endl;
            if (j==-1 ) {//scattered out of sample
                //    cout<<i<<' '<<p0.o.theta<<endl;
                // cout<<i<<' '<<tP0.theta()<<' '<<p0.o.st<<' '<<p0.o.sp<<endl;
                if (p0.scattered && p0.o.sp>=0. && fabs(p0.o.get_theta()-0.3*M_PI/180.)<s5h) dist_theta.at((unsigned int) (p0.o.sp*itheta_steps)) += pow(compton_ratio_factor,(int) p0.scattered);
            }
        }
        getrusage(RUSAGE_SELF, &r_end); //get running time
        l3+=l2;
        cout<<ii<<' '<<dist_theta.at((dist_theta.size()-1)>>1)<<' '<<l3<<' '<<l2/(r_end.ru_utime.tv_sec -r_start.ru_utime.tv_sec+ double(1e-6)*(r_end.ru_utime.tv_usec -r_start.ru_utime.tv_usec ) )<<" P/s\n";
        r_start=r_end;
        /*
        */
        out1.open(fn.c_str());
        double fac2=1./(s5h*pixelHeight*l3);
        double dx=1./itheta_steps;
        double x0=0.5*dx;
        for (unsigned int i2=0;x0<1.;i2++){
            out1<<x0<<' '<<dist_theta.at(i2)*fac2<<' '<<dist_theta.at(i2)<<' '<<l3<<endl;
            x0+=dx;
        }
        out1<<dist_theta.size()/itheta_steps<<' '<<dist_theta.at((dist_theta.size()-1)>>1)*fac2<<' '<<s5h<<' '<<pixelHeight<<endl;
        out1.close();
        ii++;
    }





    return(0);
}
