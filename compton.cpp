#include "compton.h"
using namespace std;

double e_new(double theta0,void *params) // compton formula
{
    double e0= * (double*) params;
    return(1./(1.+e0*(1-cos(theta0))));
}

double dSigma_dOmega(double theta0,void *params) // Klein-Nishina formula, 1/2 r_e^2 excluded
{
    double epsilon1=e_new(theta0,params);
    double sin1=sin(theta0);
    double sin2=sin1*sin1;
    return( sin1*epsilon1*(1. - epsilon1*(sin2-epsilon1)));
}


//double thetaDistribution::dummy(void* a, double x){return (static_cast<thetaDistribution*>(a)->p0(double x));}

thetaDistribution::thetaDistribution(double e_photon) //photon energy
{
    e0=e_photon/E_electron; // accept input in KeV
    int w_size=8096;
    gsl_integration_workspace * w0 = gsl_integration_workspace_alloc (w_size);
    gsl_function F;
    F.function = &dSigma_dOmega;
    double alpha;
    F.params = (void*) (&e0);
    double err0=1e-7;
    //gsl_integration_qags (&F, 0, M_PI, err0, err0, w_size, w0, &ans, &alpha);
    //sigma_total= 2*M_PI*ans;
    double x=0,y=0,dy,dx=M_PI/int_steps;
    thetai.resize(int_steps+1);
    fi.resize(int_steps+1);
    thetai.at(0)=x;
    fi.at(0)=y;
    double ans=0.;

//#pragma omp parallel for private(i,x,F,w0,err0,dy,alpha) schedule(static)
    for (int i=1;i<=int_steps;i++) {
        x=i*dx;
        gsl_integration_qags (&F, x-dx, x, err0, err0, w_size, w0, &dy, &alpha);
        thetai.at(i)=x;
        fi.at(i)=dy;
        ans += dy;
    }
    sigma_total= 2*M_PI*ans;
    ans=1./ans;
    y=0.;
    for (int i=1;i<=int_steps;i++) {
        y += fi.at(i);
        fi.at(i)=y*ans;
        //    cout<<thetai.at(i)<<' '<<fi.at(i)<<endl;
    }
    gsl_integration_workspace_free(w0);//free integration work space
    /* interpolation */
    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_cspline, fi.size());
    double *pxi=&(fi[0]),*pyi=&(thetai[0]);
    gsl_spline_init (spline, pxi, pyi, fi.size());
}

double thetaDistribution::theta(void) //generate polar angle theta
{
    return( gsl_spline_eval (spline, (double) random()/(RAND_MAX+(double) 1.), acc));
}
