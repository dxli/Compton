#include<cstdlib>
int mcd(int i,int j)
//maximum common divider
{
    if (!i || !j) return(0);
    if (abs(i)==1 || abs(j)==1) return(1);
    int k,i1,j1;
    if (i<j) {
        i1=j;
        j1=i;
    }else{
        i1=i;
        j1=j;
    }
    do{
        k=i1%j1;
        i1=j1;
        j1=k;
    }while (k);
    return(i1);
}


