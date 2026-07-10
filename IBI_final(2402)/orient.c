#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
double calc_orient(double **coords, int Nbb, int nchain, int natom_perchain, double box_size, int beg);
// calculate the end-to-end distance of all backbones
double calc_orient(double **coords, int Nbb, int nchain, int natom_perchain, double box_size, int beg)
{
    double res = 0, sqrr;
    double rx1, ry1, rz1;
    double rx2, ry2, rz2;
    int i, j, bb_beg, bb_end;
    double dist, cos_theta = 0, count_theta = 0, P2, norm_i, norm_j;
    double *vector_x, *vector_y, *vector_z;
    vector_x = (double *) malloc(nchain * sizeof(double));
    vector_y = (double *) malloc(nchain * sizeof(double));
    vector_z = (double *) malloc(nchain * sizeof(double));
    for (i = 0; i < nchain; ++i)
    {
        bb_beg = beg+i * natom_perchain;
        bb_end = bb_beg + Nbb - 1;
        // printf("%d %d\n", bb_beg, bb_end);
        rx1 = coords[bb_beg][0];
        ry1 = coords[bb_beg][1];
        rz1 = coords[bb_beg][2];
        rx2 = coords[bb_end][0];
        ry2 = coords[bb_end][1];
        rz2 = coords[bb_end][2];
        vector_x[i] = rx2 - rx1;
        vector_y[i] = ry2 - ry1;
        vector_z[i] = rz2 - rz1;
        sqrr = sqrt(vector_x[i] * vector_x[i] + vector_y[i] * vector_y[i] + vector_z[i] * vector_z[i]);
        // printf("%lf %lf %lf\n", vector_x[i], vector_y[i], vector_z[i]);
        // printf("%lf\n", sqrr);
    }
    for (i = 0; i < nchain-1; ++i)
    {
        for(j = i+1; j < nchain; ++j)
        {
            norm_i = sqrt(vector_x[i] * vector_x[i] + vector_y[i] * vector_y[i] + vector_z[i] * vector_z[i]);
            norm_j = sqrt(vector_x[j] * vector_x[j] + vector_y[j] * vector_y[j] + vector_z[j] * vector_z[j]);

            if (norm_i > 1e-6 && norm_j > 1e-6) 
            {  // Avoid division by zero
                cos_theta += (vector_x[i] * vector_x[j] + vector_y[i] * vector_y[j] + vector_z[i] * vector_z[j]) / (norm_i * norm_j);
                count_theta += 1;
                printf("%lf\n",(vector_x[i] * vector_x[j] + vector_y[i] * vector_y[j] + vector_z[i] * vector_z[j]) / (norm_i * norm_j));
            }
        }
    }
    P2 = 1.5 * (double) cos_theta / count_theta - 0.5;
    free(vector_x);free(vector_y);free(vector_z);
    return P2;
}

int main (int argc, const char * argv[])
{

    int i, j, t, Nbb, Nsc, f_branch, natom_perchain, natom_total, nbonds_total, pm_type, nint, nchain, nchain1, nchain2, tmax, nbin, bin, bb_perISC, maxsite_1d, nsites;
    double Nbar, density, rho_norm, bondlen, chi, kappa, kspring, box_size, grid_size, box_volume, dr, dx, dy, dz, dxo, dyo, dzo, rx, ry, rz, rxo, ryo, rzo, comx, comy, comz, r, r_low, r_high, shell_volume, number_density, normalization, gr;    
    tmax = 50;

    int *time = (int *)calloc(tmax, sizeof(int));
    double *orient_time = (double *)calloc(tmax, sizeof(double));


    FILE *outputfile, *traj, *orient;
    outputfile = fopen(argv[1], "r");
    fscanf(outputfile, "Nbb = %d\n", &Nbb);
    fscanf(outputfile, "Nsc = %d\n", &Nsc);
    fscanf(outputfile, "f_branch = %d\n", &f_branch);
    fscanf(outputfile, "natom_perchain = %d\n", &natom_perchain);
    fscanf(outputfile, "natom_total = %d\n", &natom_total);
    fscanf(outputfile, "nbonds_total = %d\n", &nbonds_total);
    fscanf(outputfile, "pm_type = %d\n", &pm_type);
    fscanf(outputfile, "nint = %d\n", &nint);
    fscanf(outputfile, "Nbar = %lf\n", &Nbar);
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
    printf("done");
    
    double **coords = (double **)calloc(natom_total, sizeof(double *));
    for (i = 0; i < natom_total; ++i)
    {
        coords[i] = (double *)calloc(3, sizeof(double));
    }
    int *particle_type = (int *)calloc(natom_total, sizeof(int));
    nchain = natom_total / natom_perchain;

    traj = fopen(argv[2], "r"); 
    for (t=0; t<tmax; t++)
    {
        // read the coordinates of the particles
        fscanf(traj, "%d\n%d 1.500000\n", &natom_total, &time[t]);
        for (i = 0; i < natom_total; ++i)
        {
            fscanf(traj, "%d %lf %lf %lf\n", &particle_type[i], &coords[i][0], &coords[i][1], &coords[i][2]);
            // if(t%5 == 0 && i < 10) printf("%d %lf %lf %lf\n", particle_type[i], coords[i][0], coords[i][1], coords[i][2]);
        }
    }
    
    orient = fopen("orient.txt", "w");
    for(t=(tmax-10); t<tmax; t++)
    {
        printf("t = %d\n", t);
        orient_time[t] = calc_orient(coords, Nbb, nchain, natom_perchain, box_size, 0);
        fprintf(orient, "%d %lf\n", time[t], orient_time[t]);
    }
    fclose(orient);


    free(time);free(orient_time);free(particle_type);
    for (i = 0; i < natom_total; ++i) free(coords[i]);
    free(coords);
    return 0;
}



