#include "header.h"

long initRan()
{
   //This will hopefully allow us to have a unique seed even if executed multiple times a second-Got from Mike
   //http://stackoverflow.com/questions/322938/recommended-way-to-initialize-srand
   unsigned long a = clock();
   unsigned long b = time(NULL);
   unsigned long c = getpid();
   a=a-b;  a=a-c;  a=a^(c >> 13);
   b=b-c;  b=b-a;  b=b^(a << 8);
   c=c-a;  c=c-b;  c=c^(b >> 13);
   a=a-b;  a=a-c;  a=a^(c >> 12);
   b=b-c;  b=b-a;  b=b^(a << 16);
   c=c-a;  c=c-b;  c=c^(b >> 5);
   a=a-b;  a=a-c;  a=a^(c >> 3);
   b=b-c;  b=b-a;  b=b^(a << 10);
   c=c-a;  c=c-b;  c=c^(b >> 15);
   return c%1000000000; //careful here.  Another 0 might break the ran1 (long long instead of just long)
}

float ran1(long *idum)
{
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    
    if(*idum <= 0)
    {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2 = (*idum);
        for(j=NTAB+7;j>=0;--j)
        {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if(*idum<0) *idum+=IM1;
            if(j<NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if(*idum<0) *idum += IM1;
    k=idum2/IQ2;
    if(*idum<0) idum2+= IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if(iy<1) iy += IMM1;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

float gasdev(long *idum)
{
    float ran1(long *idum);
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;
    
    if (*idum < 0) iset = 0;
    if (iset == 0)
    {
        do
        {
            v1 = 2.0*ran1(idum)-1.0;
            v2 = 2.0*ran1(idum)-1.0;
            rsq = v1*v1+v2*v2;
        }
        while(rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    }
    else
    {
        iset = 0;
        return gset;
    }
}
