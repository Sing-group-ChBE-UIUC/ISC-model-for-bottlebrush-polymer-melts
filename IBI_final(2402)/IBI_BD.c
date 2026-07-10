#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

pid_t getpid(void);
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NDIV (1+IMM1/NTAB)
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define PI 3.14159265359

void cholesky (double *mm, int Chainlength, double *cc, int t);
float ran1(long *idum);
float gasdev(long *idum);
long initRan();
float stochastic();

int main (int argc, const char * argv[])
{
	//VARIABLES!
    int cfg, Nbb, Nsc, fbranch, Ntot, nchain, Nisc, Natmtype, Nbonds, Nbndtype, Nangles, Nangletype, index, type, molecule;
    int dint1, dint2, dint3, low_index, high_index, Ntot_gr;
    int i, j, k, t, tmax, teq, maxiter, dmax, iteration, runtype, Nsave, thermosave, Nbin, bin, thermocount;
    double dt, p, rcut, kappa, kappa_bend, Lx, Ly, Lz, alpha, merit_cutoff, binsize;
    double dx, dy, dz, dr, dx1, dy1, dz1, dx2, dy2, dz2, coeff, cos_theta, r, r1, r2, rr, rr1, rr2, Fs, Fb; 
    double U, Utot, F, Ftot, Unb, Unbtot, bondlenavg, bondlen;
    double r_low, r_high, shell_volume, number_density, norm, merit, merit_prev,target_integral, deviation_integral;
    thermocount = 0;
    Utot = 0.0;
    Ftot = 0.0;
    Unbtot = 0.0;
    bondlenavg = 0.0;
    iteration = 0;
    U = 0.0;
    F = 0.0;
    Unb = 0.0;
    bondlen = 0.0;

	//File Input
	FILE *inputfile, *cfg_ini, *pot_ini, *gr_ini;
	inputfile = fopen("Input.txt", "r");
    cfg_ini = fopen("final0.cfg", "r");
    pot_ini = fopen("pot0.table", "r");
    gr_ini = fopen("RDF0.txt", "r");
	time_t theTime;
	time(&theTime);  
    //Read Input File  
    fscanf(inputfile, "cfg_yes = %d\n", &cfg);
    fscanf(inputfile, "Nbb = %d\n", &Nbb);
    fscanf(inputfile, "Nsc = %d\n", &Nsc);
    fscanf(inputfile, "fbranch = %d\n", &fbranch);
    fscanf(inputfile, "Nisc = %d\n", &Nisc);
    fscanf(inputfile, "kappa = %lf\n", &kappa);
    fscanf(inputfile, "kappa_bend = %lf\n", &kappa_bend);
    fscanf(inputfile, "dt = %lf\n", &dt);
    fscanf(inputfile, "rcut = %lf\n", &rcut);
    fscanf(inputfile, "tmax = %d\n", &tmax);
    fscanf(inputfile, "teq = %d\n", &teq);
    fscanf(inputfile, "maxiter = %d\n", &maxiter);
    fscanf(inputfile, "Nsave = %d\n", &Nsave);
    fscanf(inputfile, "thermosave = %d\n", &thermosave);
    fscanf(inputfile, "Nbin = %d\n", &Nbin);
    fscanf(inputfile, "alpha = %lf\n", &alpha);
    fscanf(inputfile, "merit_cutoff = %lf\n", &merit_cutoff);
    fclose(inputfile);
    dmax = Nbin;
    printf("%lf", rcut);

	//Read cfg file
    if(cfg ==1) 
	{
		fscanf(cfg_ini, "#Bottlebrush coarse-grained configuration\n");
        fscanf(cfg_ini, "%d atoms\n%d atom types\n%d bonds\n%d bond types\n%d angles\n%d angle types\n", &Ntot, &Natmtype, &Nbonds, &Nbndtype, &Nangles, &Nangletype);
        fscanf(cfg_ini, "0.0 %lf xlo xhi\n0.0 %lf ylo yhi\n0.0 %lf zlo zhi\n0 0 0 xy xz yz\n\n", &Lx, &Ly, &Lz);
        fscanf(cfg_ini, "Atoms # angle\n\n");
    }
    else if(cfg ==0)
    {}
    printf("structure read");
    nchain = (int)Ntot/Nisc;

    //File output
    FILE *BDout, *stats, *gri_IBI, *pmfi_IBI;
    char* str1 = malloc(sizeof(char)*30);
    char* str2 = malloc(sizeof(char)*30);
    char* str3 = malloc(sizeof(char)*30);
    char* str4 = malloc(sizeof(char)*30);
    sprintf(str1, "BD_IBI.txt");
    sprintf(str2, "thermo_IBI");
    sprintf(str3, "gri_IBI");
    sprintf(str4, "pmfi_IBI");

    //Memory allocation
    double *rx_ini = calloc(Ntot, sizeof(double));
    double *ry_ini = calloc(Ntot, sizeof(double));
    double *rz_ini = calloc(Ntot, sizeof(double));
    double *rx = calloc(Ntot, sizeof(double));
    double *ry = calloc(Ntot, sizeof(double));
    double *rz = calloc(Ntot, sizeof(double));
	double *fx = calloc(Ntot, sizeof(double));
    double *fy = calloc(Ntot, sizeof(double));
    double *fz = calloc(Ntot, sizeof(double));

    double *rp = calloc(Nbin, sizeof(double));
    double *pmfo = calloc(Nbin, sizeof(double));
    double *mfo = calloc(Nbin, sizeof(double));
    double *ro = calloc(Nbin, sizeof(double));
    double *gro = calloc(Nbin, sizeof(double));
    
    double *gr = calloc(Nbin, sizeof(double));
    double *gri = calloc(Nbin, sizeof(double));
    double *pmfi = calloc(Nbin, sizeof(double));
    double *pmfi_1 = calloc(Nbin, sizeof(double));
    double *mfi = calloc(Nbin, sizeof(double));
    double *mfi_1 = calloc(Nbin, sizeof(double));

    //Initialize Random Matrix
	double *RR = calloc(Ntot*3, sizeof(double));
	p = sqrt(2.0*dt);

	//Initialize Random Variable
	long *idum = malloc(sizeof(long));
	*idum = initRan();
    
    //Initialize arrays
    for(i=0; i<Nbin; ++i)
    {
        gr[i] = 0.0;
        gri[i] = 0.0;
        pmfi[i] = 0.0;
        pmfi_1[i] = 0.0;
        mfi[i] = 0.0;
        mfi_1[i] = 0.0;
    }

    //Initial structure
    if(cfg ==1) 
	{
        for(i=0; i<nchain; ++i)
		{
			for(j=0; j<Nisc; ++j)
			{
				fscanf(cfg_ini, "%d %d %d %lf %lf %lf\n", &index, &molecule, &type, &rx[i*Nisc+j], &ry[i*Nisc+j], &rz[i*Nisc+j]);
                if(molecule != i+1)
                {
                    printf("Error: Molecule number does not match chain number\n");
                    exit(1);
                }
            }
		}
        fclose(cfg_ini);
	}
	else 
	{
        // generate random structure
	}
	
    //Initial potentials 
	for(i=0; i<dmax; i++)
	{
        fscanf(pot_ini, "%d %lf %lf %lf\n", &index, &rp[i], &pmfo[i], &mfo[i]);
        if(i>0) binsize = rp[i]-rp[i-1];
    }
    fclose(pot_ini);
    printf("pot read");
    
    //Initial pmfi
    for(i=0; i<Nbin; i++)
    {
        pmfi[i] = pmfo[i];
        mfi[i] = mfo[i];
    }
	//Target g(r)
	for(i=0; i<Nbin; i++)
	{
        fscanf(gr_ini, "%lf %lf\n", &ro[i], &gro[i]);
        
        if(ro[i] != rp[i])
        {
            printf("Error: r values of gr and pmf do not match\n");
            exit(1);
        }
	}
    fclose(gr_ini);
    //Initial coordinates
    for(i=0; i<Ntot; ++i)
    {
        rx_ini[i] = rx[i];
        ry_ini[i] = ry[i];
        rz_ini[i] = rz[i];
    }
    
    //Main IBI method
    for(iteration=0; iteration<maxiter; iteration++)
	{
        for(i=0; i<Ntot; ++i)
        {
            rx[i] = rx_ini[i];
            ry[i] = ry_ini[i];
            rz[i] = rz_ini[i];
        }
        
        for(i=0; i<Nbin; ++i)
        {
            gr[i] = 0.0;
        }
        //BD simulation (!! Check equilibrium)
		for(t = 0; t<tmax; ++t)
		{
			for(i = 0; i<Ntot; ++i)
			{
				fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0;
			}
        #include "Interactions_IBI.h"
			//Set random velocities
			for(i = 0; i<Ntot*3; ++i)
			{
				RR[i] = gasdev(idum);
			}
			//Update positions based on forces, flows, random displacements
			for(i = 0; i<Ntot; ++i)
			{
                rx[i] += dt*fx[i]+p*RR[3*i];
                ry[i] += dt*fy[i]+p*RR[3*i+1];
                rz[i] += dt*fz[i]+p*RR[3*i+2];
                if(rx[i] > Lx) rx[i]-=Lx;
                if(rx[i] < 0.0) rx[i]+=Lx;
                if(ry[i] > Ly) ry[i]-=Ly;
                if(ry[i] < 0.0) ry[i]+=Ly;
                if(rz[i] > Lz) rz[i]-=Lz;
                if(rz[i] < 0.0) rz[i]+=Lz;
			}
			//Every so often, calculate g(r)
			if(t > (tmax-tmax*0.9))
			{
                for(i = 0; i<Ntot; ++i)
				{
					for(j = 0; j<Ntot; ++j)
					{
                        if(i!=j) 
                        {   
                            if((int)i/Nisc == (int)j/Nisc && abs((int)i%Nisc-(int)j%Nisc) == 1)
                            {}
                            else
                            {
                                dx = rx[i]-rx[j];
                                dy = ry[i]-ry[j];
                                dz = rz[i]-rz[j];
                                if(dx > Lx/2.0) dx-=Lx;if(dx < -Lx/2.0) dx+=Lx;
                                if(dy > Ly/2.0) dy-=Ly;if(dy < -Ly/2.0) dy+=Ly;
                                if(dz > Lz/2.0) dz-=Lz;if(dz < -Lz/2.0) dz+=Lz;
                                dr=sqrt(dx*dx+dy*dy+dz*dz);
                                bin = (int)(dr/binsize);
                                if(bin < Nbin)
                                {
                                    gr[bin]++;
                                }
                            }
                        }
					}
				}
			}	
            //Save coordinates
            if(t%Nsave == 0)
            {
                BDout = fopen(str1, "a");
                fprintf(BDout, "%d\n%d %d\n", Ntot, (teq+t), iteration);
                for(i=0; i<Ntot; ++i)
                {
                    fprintf(BDout, "%d %d %lf %lf %lf\n", i+1, 1, rx[i], ry[i], rz[i]);
                }
                fclose(BDout);
            }
            //Save thermo
            if(t%thermosave == 0)
            {
                Utot += (double)U/thermosave;
                Ftot += (double)F/thermosave;
                Unbtot += (double)Unb/thermosave;
                bondlenavg += (double)bondlen/thermosave;
                thermocount++;
                stats = fopen(str2, "a");
                fprintf(stats, "%d %lf %lf %lf %lf %lf %d\n", t, Utot/(double)thermocount, Ftot/(double)thermocount, Unbtot/(double)thermocount, bondlenavg/(double)thermocount, merit_prev, iteration);
                fclose(stats);
                U=0.0; F=0.0; Unb=0.0; bondlenavg=0.0;
            }		
        }

		///Calculate and evaluate g(r)
        gri_IBI=fopen(str3, "a");
        fprintf(gri_IBI, "Iteration %d\n", iteration);
		for(i = 0; i<Nbin; ++i)
		{
            r_low = i*binsize;
			r_high = (i+1)*binsize;
			shell_volume = 4*M_PI*r_low*r_low*binsize;
			number_density = Ntot_gr/(Lx*Ly*Lz);
			norm = (double) Ntot_gr*shell_volume*number_density*(tmax*0.9-1);
			gri[i] = (double) gr[i]/norm;
            fprintf(gri_IBI, "%lf %lf\n", ro[i], gri[i]);
		}
        fclose(gri_IBI); 
        //Merit function
        target_integral = 0.0;
        deviation_integral = 0.0;
        for(i=1; i<Nbin-1; i++)
        {
            
            if(isnan(gro[i]) || isinf(gro[i])) continue;
            else
            {
                target_integral += 0.5*(gro[i]*gro[i]+gro[i+1]*gro[i+1])*binsize;
                deviation_integral += 0.5*((gro[i]-gri[i])*(gro[i]-gri[i])+((gro[i+1]-gri[i+1])*(gro[i+1]-gri[i+1])))*binsize;
            }
        }
        merit = (double)deviation_integral/target_integral;
        printf("merit_prev %lf\n", merit_prev);
        printf("merit %lf\n", merit);

		if(merit < merit_cutoff)
		{
			printf("within cutoff, break");
            break;
		}
        else if (merit > merit_prev)
        {
            merit_prev = merit;
            printf("over damp, go inverse iteration");
        }
        else if (merit < merit_prev) 
        {
            merit_prev = merit;
            printf("under damp, go iteration");
        }
        /// IBI by comparing gtarget 
        for(i=0; i<Nbin; i++)
        {
            if(gri[i]>0 && gro[i]>0 && !isinf(gri[i]))
            {
                pmfi_1[i] = pmfi[i]-alpha*log((double)gro[i]/gri[i]);
            }
            else
            {
                pmfi_1[i] = pmfi[i];
            }
            pmfi[i] = pmfi_1[i]; //update pmf
        }
        for(i=1; i<Nbin-1; i++)
        {
            mfi_1[i] = -(double)(pmfi[i+1]-pmfi[i-1])/(2*binsize);
            mfi[i] = mfi_1[i]; //update mf
        }
        mfi_1[Nbin-1] = 0.0;
        mfi_1[0] = -(double)(pmfi[1]-pmfi[0])/(binsize);
        mfi[0] = mfi_1[0];
        mfi[Nbin-1] = 0.0;
        
        //Save potential
        pmfi_IBI=fopen(str4, "a");
        fprintf(pmfi_IBI, "Iteration %d\n", iteration);
        for(i=0; i<Nbin; ++i)
        {
            fprintf(pmfi_IBI, "%lf %lf %lf\n", ro[i], pmfi[i], mfi[i]);
        }
        fclose(pmfi_IBI);
	}
    free(rx); free(ry); free(rz);
    free(fx); free(fy); free(fz);
    free(rp); free(pmfo); free(mfo); free(ro); free(gro);
    free(gr); free(gri); free(pmfi); free(pmfi_1); free(mfi); free(mfi_1);
    free(RR); free(idum);
    free(str1); free(str2); free(str3); free(str4);
	return 0;

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

long initRan() {
    time_t seconds;
    time(&seconds); //switched from time due to crashes... unsure why...
    return -1*(unsigned long)(seconds); //Dividing by 100 keeps within correct bounds? not sure why this works
}







