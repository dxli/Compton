//read matrix binary data in gnuplot format
#include<iostream>
#include<fstream>
#include<vector>
using namespace std;

int main(int argc, char *argv[])
{
        if(argc<2) {
                cout<<"usage:: <matrix_file>\n";
                return(0);
        }
        ifstream in1(argv[1],ifstream::binary);
        if(! in1.is_open()){
                cout<<"usage:: <matrix_file>\n";
                return(0);
        }
        float n1;
        in1.read((char*)(&n1),sizeof(float));
        vector<float> x,y;
        x.resize( (unsigned int) ( n1 -0.5));
        y.resize( (unsigned int) ( n1 +0.5));
        cout<<x.size()<<endl;
        if(in1.fail()){
                cerr<<"incomplete data file\n";
                return 0;
        }
        in1.read((char *) (& x[0]),x.size()*sizeof(float));
        do{
                in1.read((char*) (&y[0]),sizeof(float)*y.size());
                if(in1.fail()) break;
                for(unsigned int i=0;i<x.size();i++){
                        cout<<x.at(i)<<' '<<y.at(0)<<' '<<y.at(i+1)<<endl;
                }
        }while(! in1.eof());
        return 0;
}
