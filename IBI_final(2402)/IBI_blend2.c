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
    int ntot1, ntot2, cfg, Nbb, Nsc, fbranch, Ntot, nchain, Nisc, Natmtype, Nbonds, Nbndtype, Nangles, Nangletype, index, molecule;
    int dint1, dint2, dint3, low_index, high_index, Ntot_gr;
    int i, j, k, t, tmax, teq, maxiter, dmax, iteration, runtype, Nsave, thermosave, Nbin, bin, grcount, thermocount;
    double dt, p, rcut, kappa, kappa_bend, Lx, Ly, Lz, alpha, merit_cutoff, binsize;
    double dx, dy, dz, dr, dx1, dy1, dz1, dx2, dy2, dz2, coeff, cos_theta, r, r1, r2, rr, rr1, rr2, Fs, Fb; 
    double U, Utot, F, Ftot, Unb, Unbtot, bondlenavg, bondlen;
    double r_low, r_high, shell_volume, number_density, norm, merit, merit_prev,target_integral, deviation_integral;
    grcount = 0;
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
	FILE *inputfile, *cfg_ini, *pot_ini, *potab_ini, *gro_ini, *grab_ini;
	inputfile = fopen("Input.txt", "r");
    cfg_ini = fopen("final0.cfg", "r");
    pot_ini = fopen("pot0.table", "r");
    gro_ini = fopen("RDFAA.txt", "r");
    grab_ini = fopen("RDFAB.txt", "r");
    potab_ini = fopen("potab0.table", "r");
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
    ntot1= 7*56;
    ntot2 = 7*152;
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
    FILE *BDout, *stats, *gri_IBI_0, *gri_IBI_AB, *pmfi_IBI_AB;
    char* str1 = malloc(sizeof(char)*30);
    char* str2 = malloc(sizeof(char)*30);
    char* str3 = malloc(sizeof(char)*30);
    char* str4 = malloc(sizeof(char)*30);
    char* str5 = malloc(sizeof(char)*30);
    sprintf(str1, "BD_IBI.txt");
    sprintf(str2, "thermo_IBI.txt");
    sprintf(str3, "gri_IBI_0.txt");
    sprintf(str4, "gri_IBI_AB.txt");
    sprintf(str5, "pmfi_IBI_AB.txt");


    //Memory allocation
    double *rx = calloc(Ntot, sizeof(double));
    double *ry = calloc(Ntot, sizeof(double));
    double *rz = calloc(Ntot, sizeof(double));
	double *fx = calloc(Ntot, sizeof(double));
    double *fy = calloc(Ntot, sizeof(double));
    double *fz = calloc(Ntot, sizeof(double));
    int *type = calloc(Ntot, sizeof(int));

    double *rp = calloc(Nbin, sizeof(double));
    double *pmfo = calloc(Nbin, sizeof(double));
    double *pmfo_bb = calloc(Nbin, sizeof(double));
    double *pmfo_ab = calloc(Nbin, sizeof(double));
    double *mfo = calloc(Nbin, sizeof(double));
    double *mfo_bb = calloc(Nbin, sizeof(double));
    double *mfo_ab = calloc(Nbin, sizeof(double));
    double *ro_AA = calloc(Nbin, sizeof(double));
    double *ro_AB = calloc(Nbin, sizeof(double));
    double *gro_AA = calloc(Nbin, sizeof(double));
    double *gro_AB = calloc(Nbin, sizeof(double));
    
    double *gr_AA = calloc(Nbin, sizeof(double));
    double *gr_AB = calloc(Nbin, sizeof(double));
    double *gri_AA = calloc(Nbin, sizeof(double));
    double *gri_AB = calloc(Nbin, sizeof(double));
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
        gr_AA[i] = 0.0;
        gr_AB[i] = 0.0;
        gri_AA[i] = 0.0;
        gri_AB[i] = 0.0;
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
				fscanf(cfg_ini, "%d %d %d %lf %lf %lf\n", &index, &molecule, &type[i*Nisc+j], &rx[i*Nisc+j], &ry[i*Nisc+j], &rz[i*Nisc+j]);
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
        fscanf(potab_ini, "%d %lf %lf %lf\n", &index, &rp[i], &pmfo_ab[i], &mfo_ab[i]);
        pmfo_bb[i] = pmfo_ab[i]; mfo_bb[i] = mfo_ab[i];
        pmfo_ab[i] = (pmfo[i]+pmfo_bb[i])/2.0; mfo_ab[i] = (mfo[i]+mfo_bb[i])/2.0;
        if(i>0) binsize = rp[i]-rp[i-1];
    }
    fclose(pot_ini);
    fclose(potab_ini);
    printf("pot read");

	//Target g(r)
	for(i=0; i<Nbin; i++)
	{
        fscanf(gro_ini, "%lf %lf\n", &ro_AA[i], &gro_AA[i]);
        
        if(ro_AA[i] != rp[i])
        {
            printf("%lf %lf\n", ro_AA[i], rp[i]);
            printf("Error: r values of gr and pmf do not match\n");
            exit(1);
        }
	}
    fclose(gro_ini);
    for(i=0; i<Nbin; i++)
    {
        fscanf(grab_ini, "%lf %lf\n", &ro_AB[i], &gro_AB[i]);
        
        if(ro_AB[i] != rp[i])
        {
            printf("Error: r values of gr and pmf do not match\n");
            exit(1);
        }
    }

	////Simulation loop - Equilibrium+ final 1000 steps get g(r)
	for(t = 0; t<teq; ++t)
	{
		for(i = 0; i<Ntot; ++i)
		{
			fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0;
		}

	#include "Interactions_blend2.h"

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
		//Fianl 1000 steps, calculate g(r)
		if(t>(teq-(int)(teq*0.3)))
		{
            for(i = 0; i<Ntot; ++i)
			{
				for(j = 0; j<Ntot; ++j)
				{
					if(i!=j && type[i]==1 && type[j]==1) 
					{   
                        if((int)i/Nisc == (int)j/Nisc && abs((int)i%Nisc-(int)j%Nisc) == 1)
                        {}
                        else
                        //if((int)((int)i%Nisc)%2==0 && (int)((int)j%Nisc)%2==0)
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
                                gr_AA[bin]++;
                                //grcount++;
                            }
                        }
					}
				}
			}

            for(i = 0; i<Ntot; ++i)
            {
                for(j = 0; j<Ntot; ++j)
                {
                    if(type[i]!=type[j]) 
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
                            gr_AB[bin]++;
                            //grcount++;
                        }
                    }
                }
            }
            
            //grcount++;
            //printf("grcount = %d\n", grcount);
		}
        //Save coordinates
        if(t%Nsave == 0)
        {
            BDout=fopen(str1, "a");
            fprintf(BDout, "%d\n%d %d\n", Ntot, t, iteration);
            for(i=0; i<Ntot; ++i)
            {
                fprintf(BDout, "%d %d %d %lf %lf %lf\n", i+1, 1, type[i], rx[i], ry[i], rz[i]);
            }
            fclose(BDout);
        }
        //Save thermo
        if(t%thermosave == 0)
        {
            Utot += U/(double)thermosave;
            Ftot += F/(double)thermosave;
            //Unbtot += Unb/(double)thermosave;
            bondlenavg += bondlen/(double)thermosave;
            thermocount++;
            stats = fopen(str2, "a");
            fprintf(stats, "%d %lf %lf %lf %d\n", t, Utot/(double)thermocount, Ftot/(double)thermocount, bondlenavg/(double)thermocount, iteration);
            fclose(stats);
            U=0.0; F=0.0; Unb=0.0; bondlen=0.0;
            // if(Utot/(double)thermocount < -1.0)
            // {
            //     printf("Error: Energy too low\n");
            //     exit(1);
            // }
        }
    }
    thermocount = 0;

	///Calculate g(r)
    //AA pair
	gri_IBI_0=fopen(str3, "a");
    fprintf(gri_IBI_0, "Iteration %d\n", iteration);
    for(i = 0; i<Nbin; ++i)
	{
		//Ntot_gr=nchain*((int)(Nisc+1)/2);
        Ntot_gr=ntot1;
        r_low = i*binsize;
		r_high = (i+1)*binsize;
		shell_volume = 4*M_PI*r_low*r_low*binsize;
		number_density = Ntot_gr/(Lx*Ly*Lz);
		norm = (double) Ntot_gr*shell_volume*number_density*(teq*0.3-1);
		gri_AA[i] = (double) gr_AA[i]/norm;
        fprintf(gri_IBI_0, "%lf %lf\n", ro_AA[i], gri_AA[i]);
	}
    fclose(gri_IBI_0);
    //AB pair
    gri_IBI_AB=fopen(str4, "a");
    fprintf(gri_IBI_AB, "Iteration %d\n", iteration);
    for(i = 0; i<Nbin; ++i)
    {
        //Ntot_gr=nchain*((int)(Nisc+1)/2);
        Ntot_gr=ntot1;
        r_low = i*binsize;
        r_high = (i+1)*binsize;
        shell_volume = 4*M_PI*r_low*r_low*binsize;
        number_density = ntot2/(Lx*Ly*Lz);
        norm = (double) Ntot_gr*shell_volume*number_density*(teq*0.3-1);
        gri_AB[i] = (double) gr_AB[i]/norm;
        fprintf(gri_IBI_AB, "%lf %lf\n", ro_AB[i], gri_AB[i]);
    }
    fclose(gri_IBI_AB);
    
	// IBI by comparing gtarget 
	for(i=0; i<Nbin; i++)
	{
		if(gro_AB[i]>0 && gri_AB[i]>0 && !isinf(gri_AB[i]))
		{
			pmfi_1[i] = pmfo_ab[i]-alpha*log((double)gro_AB[i]/gri_AB[i]);
            //printf("there is gro, gri values to compare\n");
		}
		else
		{
			pmfi_1[i] = pmfo_ab[i];
            //printf("bring pmfo value\n");
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
    pmfi_IBI_AB=fopen(str5, "a");
    fprintf(pmfi_IBI_AB, "Iteration %d\n", iteration);
    for(i=0; i<Nbin; ++i)
    {
        fprintf(pmfi_IBI_AB, "%lf %lf %lf\n", ro_AB[i], pmfi[i], mfi[i]);
    }
    fclose(pmfi_IBI_AB);

	//Main IBI method
	runtype = 0;
    //Initialize arrays
    for(iteration=0; iteration<maxiter; iteration++)
	{
        for(i=0; i<Nbin; ++i)
        {
            gr_AA[i] = 0.0;
            gr_AB[i] = 0.0;
        }
        //BD simulation (!! Check equilibrium)
		for(t = 0; t<tmax; ++t)
		{
			for(i = 0; i<Ntot; ++i)
			{
				fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0;
			}
        #include "Interactions_blend2.h"
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
                        if(i!=j && type[i]==1 && type[j]==1) 
                        {   
                            if((int)i/Nisc == (int)j/Nisc && abs((int)i%Nisc-(int)j%Nisc) == 1)
                            {}
                            else
                            //if((int)((int)i%Nisc)%2==0 && (int)((int)j%Nisc)%2==0)
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
                                    gr_AA[bin]++;
                                    //grcount++;
                                }
                            }
                        }
                    }
                }

                for(i = 0; i<Ntot; ++i)
                {
                    for(j = 0; j<Ntot; ++j)
                    {
                        if(type[i]!=type[j]) 
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
                                gr_AB[bin]++;
                                //grcount++;
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
                    fprintf(BDout, "%d %d %d %lf %lf %lf\n", i+1, 1, type[i], rx[i], ry[i], rz[i]);
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
                // if(Utot/(double)thermocount < -1.0)
                // {
                //     printf("Error: Energy too low\n");
                //     exit(1);
                // }
            }		
        }
        //thermocount = 0;
		///Calculate and evaluate g(r)
        //AA pair
        gri_IBI_0=fopen(str3, "a");
        fprintf(gri_IBI_0, "Iteration %d\n", iteration);
		for(i = 0; i<Nbin; ++i)
		{
			Ntot_gr=ntot1;
            r_low = i*binsize;
			r_high = (i+1)*binsize;
			shell_volume = 4*M_PI*r_low*r_low*binsize;
			number_density = Ntot_gr/(Lx*Ly*Lz);
			norm = (double) Ntot_gr*shell_volume*number_density*(tmax*0.9-1);
			gri_AA[i] = (double) gr_AA[i]/norm;
            fprintf(gri_IBI_0, "%lf %lf\n", ro_AA[i], gri_AA[i]);
		}
        fclose(gri_IBI_0); 
        //AB pair
        gri_IBI_AB=fopen(str4, "a");
        fprintf(gri_IBI_AB, "Iteration %d\n", iteration);
        for(i = 0; i<Nbin; ++i)
        {
            Ntot_gr=ntot1;
            r_low = i*binsize;
            r_high = (i+1)*binsize;
            shell_volume = 4*M_PI*r_low*r_low*binsize;
            number_density = ntot2/(Lx*Ly*Lz);
            norm = (double) Ntot_gr*shell_volume*number_density*(tmax*0.9-1);
            gri_AB[i] = (double) gr_AB[i]/norm/2;
            fprintf(gri_IBI_AB, "%lf %lf\n", ro_AB[i], gri_AB[i]);
        }
        fclose(gri_IBI_AB);

        //Merit function
        target_integral = 0.0;
        deviation_integral = 0.0;
        for(i=1; i<Nbin-1; i++)
        {
            if(isnan(gro_AB[i]) || isinf(gro_AB[i])) continue;
            else
            {
                target_integral += 0.5*(gro_AB[i]*gro_AB[i]+gro_AB[i+1]*gro_AB[i+1])*binsize;
                deviation_integral += 0.5*((gro_AB[i]-gri_AB[i])*(gro_AB[i]-gri_AB[i])+((gro_AB[i+1]-gri_AB[i+1])*(gro_AB[i+1]-gri_AB[i+1])))*binsize;
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
            runtype = 0;
            merit_prev = merit;
            printf("over damp, go iteration, quit if it last long");
        }
        else if (merit < merit_prev) 
        {
            runtype = 0;
            merit_prev = merit;
            printf("under damp, go iteration");
        }
        /// IBI by comparing gtarget 
        // if(runtype == 0)
        // {
        for(i=0; i<Nbin; i++)
        {
            if(gri_AB[i]>0 && gro_AB[i]>0 && !isinf(gri_AB[i]))
            {
                pmfi_1[i] = pmfi[i]-alpha*log((double)gro_AB[i]/gri_AB[i]);
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
        mfi_1[0] = mfi_1[1]+(rp[0]-rp[1])*(mfi_1[2]-mfi_1[1])/(double)(rp[2]-rp[1]);
        mfi[0] = mfi_1[0];
        mfi[Nbin-1] = 0.0;
        // }
        // else if (runtype == 1)
        // {}
        //Save potential
        pmfi_IBI_AB=fopen(str5, "a");
        fprintf(pmfi_IBI_AB, "Iteration %d\n", iteration);
        for(i=0; i<Nbin; ++i)
        {
            fprintf(pmfi_IBI_AB, "%lf %lf %lf\n", ro_AB[i], pmfi[i], mfi[i]);
        }
        fclose(pmfi_IBI_AB);
	}
    //Memory deallocation
    free(str1); free(str2); free(str3); free(str4); free(str5);
    free(rx); free(ry); free(rz); free(fx); free(fy); free(fz); free(type);
    free(rp); free(pmfo); free(pmfo_ab); free(mfo); free(mfo_ab); free(ro_AA); free(ro_AB); free(gro_AA); free(gro_AB);
    free(gr_AA); free(gr_AB); free(gri_AA); free(gri_AB); free(pmfi); free(pmfi_1); free(mfi); free(mfi_1);
    free(RR);
    free(idum);
    return 0;

}

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







