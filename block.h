#include <cstring>
#include <gsl/gsl_matrix.h>

#define cellDiv 200
double E_to_l(double );
class EulerAngles{
public:
    double ct,st,cp,sp;
    double sp0,cp0,ct0;
    //double theta,phi;
    EulerAngles();
    EulerAngles(double);
    EulerAngles(double,double);
    void addTo(EulerAngles);
    double get_theta(),get_phi();
};
class Coordinates{
public:
    double x,y,z;
    Coordinates();
    Coordinates(double,double,double);
    Coordinates(double,double);
    Coordinates(EulerAngles);
    double normxy();
    void addTo(Coordinates);
};
Coordinates operator *(double ,Coordinates );
Coordinates operator +(Coordinates ,Coordinates );

EulerAngles operator +(EulerAngles , EulerAngles );



class photon{
private:
    double r0,cp0,sp0;
public:
    double muPE,muC,muTotal,imuTotal;
    double rstep,R2;
    int compton_ratio;
    //double travel;
    //unsigned int scattered;
    double sampleL,muPEdr,muCdr;
    unsigned int scattered;
    Coordinates r;
    EulerAngles o;
    photon();
    photon(double);
    int init();
    int propagate(double);
    int initRef(double);
    int initRefBuried(double);
    int propagateRef(double);
    int propagatePass(double);
    int initPass();
};
