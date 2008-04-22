#include"block.h"
#include<cmath>
#include<iostream>
using namespace std;

EulerAngles::EulerAngles()
{
    //do nothing
}
EulerAngles::EulerAngles(double theta0)
//fixed theta0, and free rotation around the same theta
{
    double phi0=random()*2.*M_PI/(RAND_MAX+1.0);
    double st0=sin(theta0);
    sp=st0*sin(phi0);
    cp=sqrt(1.-sp*sp);
    ct=cos(theta0)/cp;
    st=cos(phi0)*st0/cp;
    //phi=asin(sp); //no need
    //theta=acos(ct);
}

EulerAngles::EulerAngles(double theta0,double phi0)
{
    //theta=theta0;
    //phi=phi0;
    cp=cos(phi0),sp=sin(phi0);
    ct=cos(theta0),st=sin(theta0);
}
EulerAngles operator +(EulerAngles a, EulerAngles b)
//transform to the original coordinates, add the rotation
{
    EulerAngles o;
    o.sp=b.cp*b.ct*a.sp+a.cp*b.sp;
    o.cp=sqrt(1-o.sp*o.sp);
    o.ct=(a.cp*b.cp*a.ct*b.ct -a.ct*a.sp*b.sp-b.cp*a.st*b.st)/o.cp;
    o.st=(a.cp*b.cp*a.st*b.ct -a.st*a.sp*b.sp+b.cp*a.ct*b.st)/o.cp;
//o.phi=asin(o.sp); // no need to convert it back, we only do this after an photon scattered out of the sample
//o.theta=acos(o.ct);
    return(o);
}

void EulerAngles::addTo(EulerAngles b)
//transform to the original coordinates, add the rotation
//need to speed up this
{
    sp0=b.cp*b.ct*sp+cp*b.sp;
    cp0=sqrt(1-sp0*sp0);
    ct0=(cp*b.cp*ct*b.ct -ct*sp*b.sp-b.cp*st*b.st)/cp0;
    st=(cp*b.cp*st*b.ct -st*sp*b.sp+b.cp*ct*b.st)/cp0;
    ct=ct0;
    cp=cp0;
    sp=sp0;
}

double EulerAngles::get_theta(){
    return acos(ct);
}
double EulerAngles::get_phi(){
    return asin(sp);
}
Coordinates::Coordinates(){
    //do nothing
}

Coordinates::Coordinates(double x0,double y0,double z0)
{
    x=x0;
    y=y0;
    z=z0;
}

Coordinates::Coordinates(EulerAngles a)
{
    x=a.ct;//cos(a.theta);
    y=a.st;//sin(a.theta);
    z=a.cp;//cos(a.phi);
    y*=z;
    x*=z;
    z=a.sp;//sin(a.phi);

}
Coordinates::Coordinates(double theta0,double phi0)
{
    x=cos(theta0);
    y=sin(theta0);
    z=cos(phi0);
    y*=z;
    x*=z;
    z=sin(phi0);
}

double Coordinates::normxy()
{
    return(x*x+ y*y);
}
void Coordinates::addTo(Coordinates b)
{
    x += b.x;
    y += b.y;
    z += b.z;
}

Coordinates operator +(Coordinates a,Coordinates b)
{
    Coordinates o=a;
    o.x += b.x;
    o.y += b.y;
    o.z += b.z;
    return (o);
}
Coordinates operator *(double a,Coordinates b)
{
    Coordinates o=b;
    o.x = a*b.x;
    o.y = a*b.y;
    o.z = a*b.z;
    return (o);
}

ostream & operator << ( ostream& os,Coordinates b)
{
        os<<b.x<<' '<<b.y<<' '<<b.z;
        return os;
}

photon::photon()
{
    muPE= 1.16 * 2.33; // in cm^-1, PhotoElectric
    muC= 0.15 * 2.33; // in cm^-1, compton
    muTotal=1.31 * 2.33; // in cm^-1, total
   // compton_ratio=(int) ( RAND_MAX*muC/muTotal);
    imuTotal=1./muTotal;
    sampleL=-2.54*0.25;
    R2= (1+1e-8)*sampleL* sampleL;
}

photon::photon(double diameter)
{
    muPE= 1.16 * 2.33; // in cm^-1, PhotoElectric
    muC= 0.15 * 2.33; // in cm^-1, compton
    muTotal=1.31 * 2.33; // in cm^-1, total
   // compton_ratio=(int) ( RAND_MAX*muC/muTotal);
    imuTotal=1./muTotal;
    sampleL=- 2.54*diameter;
    R2= (1+1e-8)*sampleL* sampleL;
}

int photon::init()
{
        scattered=0;
    r.y= - sampleL- log(random()/(RAND_MAX+1.0)+1e-16)*imuTotal;
    if (r.y>=sampleL) return(1);
    r.x=r.z=0.;
    //r=Coordinates(0.,r0,0.);
    //o=EulerAngles(0.,0.);
    o.cp=o.ct=1.;
    o.sp=o.st=0.;
    return(0);
}
int photon::initRef(double phi0)
{
        scattered=0;
    r0= log((RAND_MAX+1.0)/(random()+(long int) 1))*imuTotal;
cp0=cos(phi0);
r.y=r0*cp0;
    if (r.y>=sampleL) return(1);
sp0=sin(phi0);
r.z=r0*sp0;
r.x=0.;
    //r=Coordinates(0.,r0*cos(phi0),r0*sin(phi0));
    o.ct=1.;o.st=0.;
    o.cp=cp0;o.sp=sp0;
    //o=EulerAngles(0.,phi0);
    return(0);
}

int photon::initPass()
{//Pass through
    scattered=0;
    r.z=0.;
    r.x=sampleL;
    r.y=0.;
    o.cp=1.;o.sp=0.;
    o.ct=1.;o.st=0.;
    //o=EulerAngles(0.,phi0);
    return(0);
}

int photon::initGixos()
{//Pass through
    scattered=0;//corritical angle 0.06903 degree, at Si/Ga interface 0.4132 angstrom
    r.z=0.06*M_PI/180*sampleL;
    r.x=sampleL;
    r.y=0.;
    o.cp=1.;o.sp=-0.00105;
    o.ct=1.;o.st=0.;
    //o=EulerAngles(0.,phi0);
    return(0);
}

int photon::propagatePass(double theta0)
        // pass through
{
    r.addTo(log((RAND_MAX+1.0)/(random()+(long int) 1))*imuTotal*Coordinates(o)); //propagate according to exponential decay
    if (r.normxy() > R2) return(-1); // scattered out of solid sample
   // if(random()>compton_ratio) return(0);// not compton scattering
            scattered++;
        o.addTo(EulerAngles(theta0));
    return(1);
}

int photon::initRefBuried(double sphi0)
{//Reflectivity, buried interface
    scattered=0;       
    o.sp=sphi0;
    o.cp=sqrt(1.-sphi0*sphi0);
    r.z=sampleL*sphi0/o.cp;
    r.x=sampleL;
    r.y=0.;
    o.ct=1.;o.st=0.;
    //o=EulerAngles(0.,phi0);
    return(0);
}

int photon::propagateGixos(double theta0)
// scattering process for reflectivity geometry
{
        double z0=r.z;
    r.addTo(log((RAND_MAX+1.0)/(random()+(long int) 1))*imuTotal*Coordinates(o)); //propagate according to exponential decay
    if (r.z<0.) {//gixos is at total reflection
            r.z = fabs(r.z);
            o.sp = fabs(o.sp); 
    }
    if (r.normxy() > R2 ) return(-1); // scattered out of solid sample
    // if(random()>compton_ratio) return(0);// not compton scattering
    scattered++;
    o.addTo(EulerAngles(theta0));
    return(1);
}


int photon::propagateRef(double theta0)
// scattering process for reflectivity geometry
{
    r.addTo(log((RAND_MAX+1.0)/(random()+(long int) 1))*imuTotal*Coordinates(o)); //propagate according to exponential decay
    if (r.normxy() > R2 || r.z<0) return(-1); // scattered out of solid sample
    // if(random()>compton_ratio) return(0);// not compton scattering
    scattered++;
    o.addTo(EulerAngles(theta0));
    return(1);
}

int photon::propagate(double theta0)
{
    r0=random()/(RAND_MAX+1.0);
    if (r0<=muPEdr) {// absorbed
        return(0);
    }
    r0 -= muPEdr;
    if (r0<=muCdr) {// Compton scattered
            scattered++;
        o.addTo(EulerAngles(theta0));
    }
    r.addTo(rstep*Coordinates(o));
    if (r.normxy() > R2) return(-1);
    return(1);
}

double E_to_l(double en0)
        //Energy to wavelength, KeV to angstrom
        {
         double  e=1.602176e-19, h=6.626069e-34,  c=2.997925e8;
          return h*c/e*1.e7/en0;
          }
        //
