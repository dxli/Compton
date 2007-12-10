#include<iostream>
#include<fstream>
#include<time.h>
#include<sys/time.h>
#include"randomInit.h"
using namespace std;
void randomInit(char * random_buf)
{
    string sr("/dev/urandom");
    struct timeval t_start;
    gettimeofday(&t_start,NULL);
    ifstream in1(sr.c_str());
    if (! in1.is_open())
    {
        cerr<<"Can not open" <<sr<<endl;
        srandom(t_start.tv_usec);
    }
    else
    {
        in1.read(random_buf,256*sizeof(char));
        in1.close();
        initstate(t_start.tv_usec,random_buf,256);
        setstate(random_buf);
    }
}
