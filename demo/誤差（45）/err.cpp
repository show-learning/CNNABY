#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;
int main(){
    FILE *o_yuv,*e_yuv;
    unsigned char *oBuf = (unsigned char *)malloc(352*288*sizeof(unsigned char));
	unsigned char *eBuf = (unsigned char *)malloc(352*288*sizeof(unsigned char));
    o_yuv=fopen("origin.yuv", "rb");
    fread(oBuf, sizeof(unsigned char), 352*288, o_yuv);
    e_yuv=fopen("error.yuv", "rb");
    fread(eBuf, sizeof(unsigned char), 352*288, e_yuv);
    double o[352*288],e[352*288];
    double *o_tmp=(double *)malloc(1*sizeof(double)),*e_tmp=(double *)malloc(1*sizeof(double));
    int over_30=0;
    int count=0;
    for(count=0;count<352*288;count++){
        o_tmp[0]= *oBuf++;
	    o[count]=double(o_tmp[0]);
        e_tmp[0]= *oBuf++;
	    e[count]=double(e_tmp[0]);
        if(fabs(o[count]-e[count])>45){
            over_30++;
            cout<<count<<", origin-error= "<<(o[count]-e[count])<<endl;
        }
    }
    cout<<over_30;
}
