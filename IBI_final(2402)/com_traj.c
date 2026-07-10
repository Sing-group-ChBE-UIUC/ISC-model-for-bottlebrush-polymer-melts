
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

double simbox_1dgetimage(double rx, double bl);

int main (int argc, const char * argv[])
{
    //Variables
    int i, j, c, k,l, isc, excluded_beads, t, tmax, Nbb, Nsc, f_branch, natom_perchain, natom_total, nbonds_total, pm_type, nint, maxsite_1d, nsites, nchain;
    int bin, bb_beg, bb_index, isc_end, bb_isc_beg, sc_isc_beg, valid_nbb, index, b_index, a_index, nbin, N_rcut, Tot_addmidISC, sc_sum;
    int ISC_segment, Tot_ISC, natom_inISCchain, count, bb_perISC, total_bbinISC, total_beadsinISC, excluded_bb;
    double Nbar, density, rho_norm, bondlen, chi, kappa, kspring, box_size, grid_size, box_volume;
    double comx, comy, comz, rx1, ry1, rz1, rx2, ry2, rz2, dx, dy, dz, dxo, dyo, dzo, r, dr, D, normx, normy, normz, norm_len, unit_normx, unit_normy, unit_normz;
    double rx, ry, rz, rxo, ryo, rzo, tot_gr;

    //Read output file
    FILE *outputfile;
    outputfile = fopen("output.txt", "r");
    fscanf(outputfile, "nchain1 = %d\n", &nchain);
    fscanf(outputfile, "nchain2 = %d\n", &nchain);
    fscanf(outputfile, "natom_perchain1 = %d\n", &natom_perchain);
    fscanf(outputfile, "natom_perchain2 = %d\n", &natom_perchain);
    fscanf(outputfile, "natom_total = %d\n", &natom_total);
    fscanf(outputfile, "nbonds_total = %d\n", &nbonds_total);
    fscanf(outputfile, "nint = %d\n", &nint);
    fscanf(outputfile, "density = %lf\n", &density);
    fscanf(outputfile, "rho_norm = %lf\n", &rho_norm);
    fscanf(outputfile, "bondlen = %lf\n", &bondlen);
    fscanf(outputfile, "chi = %lf\n", &chi);
    fscanf(outputfile, "kappa = %lf\n", &kappa);
    fscanf(outputfile, "kspring = %lf\n", &kspring);
    fscanf(outputfile, "box_size = %lf\n", &box_size);
    fscanf(outputfile, "grid_size = %lf\n", &grid_size);
    fscanf(outputfile, "maxsite_1d = %d\n", &maxsite_1d);
    fscanf(outputfile, "nsites = %d\n", &nsites);
    fclose(outputfile);
    printf("%d\n", nchain);
    box_size = maxsite_1d * grid_size;
    box_volume = box_size * box_size * box_size;
    // nchain = density * box_volume / natom_perchain;
    tmax = 60;
    Nbb = 100;
    // //Read contour file
    // FILE *contfile_A;
    // Nbb = 100;
    // contfile_A = fopen("contour_A.txt", "r");
    // int *bbindexA = calloc(Nbb, sizeof(int));
    // int *branchlist1_A = calloc(Nbb, sizeof(int));
    // for(int i = 0; i < Nbb; ++i)
    // {
    //     fscanf(contfile_A, "%d %d\n", &bbindexA[i], &branchlist1_A[i]);
    // }
    // fclose(contfile_A);

    //ISC parameters
    ISC_segment =  atoi(argv[1]);
    Tot_ISC = ISC_segment*nchain*2;
    Tot_addmidISC = Tot_ISC+Tot_ISC -1;
    //bb_perISC = (int)(Nbb/ISC_segment);
    bb_perISC = atoi(argv[2]);
    total_bbinISC = bb_perISC*ISC_segment;
    excluded_bb = Nbb - total_bbinISC;
    printf("excluded_beads = %d\n", excluded_bb);

    //Allocate memory
    double **coords = (double **)calloc(natom_total, sizeof(double *));
    for (i = 0; i < natom_total; ++i)
    {
        coords[i] = (double *)calloc(3, sizeof(double));
    }
    double **com = (double **)calloc(Tot_ISC, sizeof(double *));
    for(i=0; i<Tot_ISC; i++)
    {
        com[i] = (double *)calloc(3, sizeof(double));
    }
    // double **com_mid = (double **)calloc(nchain*(ISC_segment-1), sizeof(double *));
    // for(i=0; i<nchain*(ISC_segment-1); i++)
    // {
    //     com_mid[i] = (double *)calloc(3, sizeof(double));
    // }
    int *particle_type = calloc(natom_total, sizeof(int));
    int *isc_type = calloc(Tot_ISC, sizeof(int));
    int *time = calloc(tmax, sizeof(int));

    //Initialize isc type 
    for(i=0; i<Tot_ISC; i++)
    {
        com[i][0] = 0.0;
        com[i][1] = 0.0;
        com[i][2] = 0.0;
    }
    for(c=0; c<nchain*2; c++)
    {
        for(i=0; i<ISC_segment; i++) isc_type[c*ISC_segment+i]=i+1;
    }
    // for(i=0; i<nchain*(ISC_segment-1); i++)
    // {
    //     com_mid[i][0] = 0.0;
    //     com_mid[i][1] = 0.0;
    //     com_mid[i][2] = 0.0;
    // }
    // for(i=0; i<Tot_addmidISC; i++)
    // {
    //     isc_type[i] = i+1;
    // }

    //Open files
    FILE *traj, *isc_comtraj;
    char *filename = malloc(sizeof(char)*35);
    char *filename2 = malloc(sizeof(char)*35);
    char* str = malloc(sizeof(char)*35);

    sprintf(filename,"traj_blend.xyz");
    sprintf(filename2,"isc_comtraj.xyz");

    traj = fopen(filename, "r"); 
    isc_comtraj = fopen(filename2, "w");

    // Read the trajectory file
    for (t=0; t<tmax; t++)
    {
        // read the coordinates of the particles
        fscanf(traj, "%d\n%d 1.500000\n", &natom_total, &time[t]);
        for (i = 0; i < natom_total; ++i)
        {
            fscanf(traj, "%d %lf %lf %lf\n", &particle_type[i], &coords[i][0], &coords[i][1], &coords[i][2]);
        }
        // print head of traj files
        fprintf(isc_comtraj, "%d\n%d 1.500000\n", Tot_ISC, time[t]);
   
    // Map ISC segment 
        // Iterate over chains
        for(c=0; c<nchain*2; c++)
        {
            bb_beg = c*(natom_perchain);
            sc_sum = 0;
            for(i=0; i<ISC_segment; i++)
            {
                // Backbone beads within ISC segment
                count=0;
                dx=0; dy=0; dz=0; dxo=0; dyo=0; dzo=0;
                bb_isc_beg = bb_beg + i*bb_perISC;
                comx=0; comy=0; comz=0;

                // consider excluded backbone 
                if(excluded_bb <0 && i == (ISC_segment-1)) valid_nbb = bb_perISC + excluded_bb;
                else valid_nbb = bb_perISC;
                for(j=0; j<valid_nbb; j++)
                {
                    rxo = coords[bb_isc_beg+j][0];
                    ryo = coords[bb_isc_beg+j][1];
                    rzo = coords[bb_isc_beg+j][2];
                    rx = simbox_1dgetimage(rxo, box_size);
                    ry = simbox_1dgetimage(ryo, box_size);
                    rz = simbox_1dgetimage(rzo, box_size);
                    dx = rx-coords[bb_isc_beg][0];
                    dy = ry-coords[bb_isc_beg][1];
                    dz = rz-coords[bb_isc_beg][2];
                    if(dx > box_size/2) rx -= box_size; if(dx < -box_size/2) rx += box_size;
                    if(dy > box_size/2) ry -= box_size; if(dy < -box_size/2) ry += box_size;
                    if(dz > box_size/2) rz -= box_size; if(dz < -box_size/2) rz += box_size;
                    com[c*ISC_segment+i][0] += rx;
                    com[c*ISC_segment+i][1] += ry;
                    com[c*ISC_segment+i][2] += rz;
                    count++;
                    // sc_isc_beg = bb_beg+Nbb+sc_sum;
                    // for(l=0; l<branchlist1_A[j]; l++)
                    // {
                    //     rxo = coords[sc_isc_beg+l][0];
                    //     ryo = coords[sc_isc_beg+l][1];
                    //     rzo = coords[sc_isc_beg+l][2];
                    //     rx = simbox_1dgetimage(rxo, box_size);
                    //     ry = simbox_1dgetimage(ryo, box_size);
                    //     rz = simbox_1dgetimage(rzo, box_size);
                    //     dx = rx-coords[bb_isc_beg][0];
                    //     dy = ry-coords[bb_isc_beg][1];
                    //     dz = rz-coords[bb_isc_beg][2];
                    //     if(dx > box_size/2) rx -= box_size; if(dx < -box_size/2) rx += box_size;
                    //     if(dy > box_size/2) ry -= box_size; if(dy < -box_size/2) ry += box_size;
                    //     if(dz > box_size/2) rz -= box_size; if(dz < -box_size/2) rz += box_size;
                    //     com[c*ISC_segment+i][0] += rx;
                    //     com[c*ISC_segment+i][1] += ry;
                    //     com[c*ISC_segment+i][2] += rz;
                    //     count++;
                    // }
                    // sc_sum += branchlist1_A[j];
                }
                //ISC coordinates as COM 
                com[c*ISC_segment+i][0] /= count;
                com[c*ISC_segment+i][1] /= count;
                com[c*ISC_segment+i][2] /= count;
                fprintf(isc_comtraj, "%d %d %lf %lf %lf\n", c+1, isc_type[c*ISC_segment+i], com[c*ISC_segment+i][0], com[c*ISC_segment+i][1], com[c*ISC_segment+i][2]);
            }                   
        }
    }

    free(coords);
    for (i = 0; i < natom_total; ++i)
    {
        free(coords[i]);
    }
    free(coords);
    for(i=0; i<Tot_ISC; i++)
    {
        free(com[i]);
    }
    free(com);
    // for(i=0; i<nchain*(ISC_segment-1); i++)
    // {
    //     free(com_mid[i]);
    // }
    free(particle_type);
    free(isc_type);
    free(time);

}
double simbox_1dgetimage(double rx, double bl)
{
    double res = 0;
    res = rx - bl * round(rx / bl);
    return res;
}



