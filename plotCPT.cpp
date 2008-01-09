//plot height structure factor, and structure factor
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<iterator>
#include<vector>
#include<algorithm>
#include<stdlib.h>
#define _GNU_SOURCE
#include<stdio.h>
#include<math.h>
#include <grace_np.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>



#define y_max 2e-4

using namespace std;

#ifndef EXIT_SUCCESS
#  define EXIT_SUCCESS 0
#endif

#ifndef EXIT_FAILURE
#  define EXIT_FAILURE -1
#endif

void my_error_function(const char *msg)
{
    fprintf(stderr, "library message: \"%s\"\n", msg);
}


int main(int argc,char *argv[])
{
double a0,a1,a2,lx,ly,n0,n1,n2,dr,a,b,c,f,*g20,*g21,*p0,
	h2[2048+2],q,
	r2[2048+2],
	q2[2048+2],
	cnr[20]={2.5e+8,2.8e+8,2.5e+8,2.5e+8,1.15e+9,1.15e+9,1.15e+9,1.2e+9},
	 *g22,d,db[15],dd,qm,dq,qr,r,v0=-1.,ans,*datab1,*datab2,*p1,*p2;
int i,i1,j,number1,number2;
FILE *f1;
char *filename1,*cp1,*cp2,string2[80],fn3[50],xc;
string fn0("si-Compton-pass-[0-9]*.txt");
    cout<<"usage:: plot si-Compton-pass-[0-9]*.txt\n";
f1=popen( (string("ls ")+fn0).c_str(),"r");
size_t buflength=0;
char *linebuf=NULL;
vector<string> fns;
while( getline( &linebuf,&buflength,f1)){ // read in filenames
        if(feof(f1)) break;
        cout<<linebuf<<endl;
        fns.push_back(string(linebuf));
}
if(fns.size()<1) {
        cout<<"No file of "<<fn0<<" found\n";
        exit(0);
}
fclose(f1);
std::sort(fns.begin(),fns.end());
// plotting
a0=0.;i=0;
if (GraceOpenVA("xmgrace", 2048, "-nosafe", "-noask","-nosigcatch","-geometry","1400x980+0+0", NULL)==-1){
        fprintf(stderr, "Can't run Grace. \n");
        exit(EXIT_FAILURE);
    }

//GracePrintf("XAXES SCALE LOGARITHMIC");
//GracePrintf("yAXES SCALE LOGARITHMIC");
GracePrintf("view xmin 0.2");
//GracePrintf("world ymax 1000");
//GracePrintf("world ymin 0.001");
GracePrintf("world xmax %g",M_PI);
GracePrintf("world xmin 0.0");
//GracePrintf("yaxis tick major 10");
//GracePrintf("xaxis tick major 10");
//GracePrintf("yaxis tick minor ticks 9");
//GracePrintf("xaxis tick minor ticks 9");
for(i=0;i<fns.size();i++){
GracePrintf("s%d on",i);
GracePrintf("s%d type xy",i);
GracePrintf("s%d linewidth 2",i);
GracePrintf("s%d symbol 0",i);
GracePrintf("s%d symbol 0",i);
GracePrintf("s%d linewidth 2",i);
GracePrintf("s%d line linestyle 1",i);
//GracePrintf("s%d legend \"S(q)\"",i);

        ifstream in1(fns.at(i).c_str());
        if (! in1.is_open()) {
                cout<<"Can not open "<<fns.at(i)<<endl;
                continue;
        }
        vector<vector<double> > xe;
        string line0;
        while( ! in1.eof()){
                std::getline(in1,line0);
                istringstream iss (line0);
                vector<double> va;
                std::copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(va));
                if (va.size()==4) xe.push_back(va);
        }
        in1.close();
        if(xe.size()>2) xe.pop_back();
        for(int i1=0;i1<xe.size();i1++){
         GracePrintf("s%d point %g,%g",i,xe.at(i1).at(0),xe.at(i1).at(1));
        }
cout<<"i= "<<i<<endl;
i++;
}


GracePrintf("legend 0.75,0.85");
GracePrintf("xaxis label \"\\+\\+2\\f{Symbol}q\"");
GracePrintf("yaxis label \"\\+\\+Intensity\"");
GracePrintf("device \"EPS\" dpi 600");
GracePrintf("device \"EPS\" OP \"color,level2,bbox:tight\"");
GracePrintf("hardcopy device \"EPS\"");
string fnout("plotCPT.eps");
GracePrintf("print to \"%s\"",fnout.c_str());
printf("print to \"%s\"\n",fnout.c_str());
GracePrintf("print");
GracePrintf("redraw");

}

