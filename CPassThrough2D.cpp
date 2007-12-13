//scattering in pass-through geometry
//pass through geometry, beam pass through Si crystal
#include "config.h"

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
    double s3h=20./640.,s3v=20./640.;// slit size
    vector<vector<float> > ScatteredI; // vector to hold 2D compton counts, single precision, for gnuplot use
    int theta_steps=(int)(M_PI/s3h+0.5);
    int phi_steps=(int)(M_PI/2/s3v+0.5);
    double itheta_steps=theta_steps/M_PI;
    double iphi_steps=phi_steps/(M_PI/2);
    ScatteredI.resize(theta_steps+2);
    for (unsigned int i=0;i<ScatteredI.size();i++) { //initialize counts to 0
            ScatteredI.at(i).resize(phi_steps+2);
        for (unsigned int j=0;j<ScatteredI.at(i).size();j++) 
            ScatteredI.at(i).at(j)=0.;
    }
    // x and y grid, for theta and phi, respectively
        vector<float> grid_x,grid_y;
        grid_x.resize(theta_steps+1);
        grid_y.resize(phi_steps);
        grix_x[0]=theta_steps;
        for(unsigned int i=1;i<grid_x.size();i++){
                grid_x.at(i)=(i+0.5)/itheta_steps;
        }
        for(unsigned int i=0;i<grid_y.size();i++){
                grid_y.at(i)=(i+0.5)/iphi_steps;
        }
        //
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
    float compton_ratio_factor=muC/muTotal; //weight factor for each scattering
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
                if (p0.scattered ) ScatteredI.at((unsigned int) (fabs(p0.o.get_phi())*iphi_steps+0.5)).at((unsigned int) (p0.o.get_theta()*itheta_steps+0.5)) += pow(compton_ratio_factor,(int) p0.scattered);
            }
        }
        getrusage(RUSAGE_SELF, &r_end); //get running time
        l3+=l2;
        cout<<ii<<' '<<ScatteredI.at(ii)<<' '<<l3<<' '<<l2/(r_end.ru_utime.tv_sec -r_start.ru_utime.tv_sec+ double(1e-6)*(r_end.ru_utime.tv_usec -r_start.ru_utime.tv_usec ) )<<" P/s\n";
        r_start=r_end;
        /*
        */
        out1.open(fn.c_str(),ofstream::binary);
        double fac2=1./(s3v*s3h*l3);
        //write matrix binary for gnuplot
        vector<float> fbuffer;
        fbuffer.resize(theta_steps);
        for (unsigned int i2=0;i2<theta_steps;i2++){ 
                if(i2) {
                        out1.write((char *) (grid_y.at(i2-1)),sizeof(float));
        vector<float>::iterator pf0=fbuffer.begin(),pd0=ScatteredI.at(i2-1).begin();
        while(pf0 != fbuffer.end()) *pf0++ = fac2* *pd0++;
                        out1.write((char *) (&fbuffer[0]),fbuffer.size()*sizeof(float));
                }else{
                        out1.write((char *) (&grid_x[0]),grid_x.size()*sizeof(float));
                }
        }
        out1.close();
        ii++;
    }
    return(0);
}
