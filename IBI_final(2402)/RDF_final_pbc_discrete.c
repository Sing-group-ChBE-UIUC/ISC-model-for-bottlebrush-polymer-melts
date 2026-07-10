
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

double simbox_1dgetimage(double rx, double bl);

int main (int argc, const char * argv[])
{
    //Variables
    int i, j, c, k,l, isc, excluded_beads, t, tmax, Nbb, Nsc, f_branch, natom_perchain, natom_total, nbonds_total, pm_type, nint, maxsite_1d, nsites, nchain;
    int bin, bb_beg, bb_index, isc_end, bb_isc_beg, sc_isc_beg, valid_nbb, index, b_index, a_index, nbin, N_rcut;
    int ISC_segment, Tot_ISC, natom_inISCchain, count, bb_perISC, total_bbinISC, total_beadsinISC, excluded_bb, shi, shj, rei, rej;
    double Nbar, density, rho_norm, bondlen, chi, kappa, kspring, box_size, grid_size, box_volume, ISC_mL, ISC_d, d;
    double rx1, ry1, rz1, rx2, ry2, rz2, dx, dy, dz, dxo, dyo, dzo, r, dr, D, normx, normy, normz, norm_len, unit_normx, unit_normy, unit_normz;
    double r_low, r_high, shell_volume, number_density, normalization, gr, comx, comy, comz;
    double rx, ry, rz, rxo, ryo, rzo, tot_gr;
    //Read output file
    FILE *outputfile;
    outputfile = fopen("output.txt", "r");
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
    box_size = maxsite_1d * grid_size;
    box_volume = box_size * box_size * box_size;
    nchain = density * box_volume / natom_perchain;
    tmax = 80;
    nbin = 596;
    dr = box_size / (2.0 * nbin);

    //ISC parameters
    ISC_segment =  atoi(argv[1]);
    Tot_ISC = ISC_segment*nchain;
    //bb_perISC = (int)(Nbb/ISC_segment);
    bb_perISC = atoi(argv[2]);
    total_bbinISC = bb_perISC*ISC_segment;
    total_beadsinISC = bb_perISC*(f_branch*Nsc+1);
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
    double **com_mid = (double **)calloc(nchain*(ISC_segment-1), sizeof(double *));
    for(i=0; i<nchain*(ISC_segment-1); i++)
    {
        com_mid[i] = (double *)calloc(3, sizeof(double));
    }
    int *particle_type = calloc(natom_total, sizeof(int));
    int *isc_type = calloc(Tot_ISC, sizeof(int));
    int *time = calloc(tmax, sizeof(int));
    int *histogram = calloc(nbin, sizeof(int));
    double *avgRDF = calloc(nbin, sizeof(double));
    double *PMF = calloc(nbin, sizeof(double));
    double *MF = calloc(nbin, sizeof(double));

    //Initialize isc type and RDF
    for(i=0; i<Tot_ISC; i++)
    {
        isc_type[i] = i;
        com[i][0] = 0.0;
        com[i][1] = 0.0;
        com[i][2] = 0.0;
    }
    for(i=0; i<nchain*(ISC_segment-1); i++)
    {
        com_mid[i][0] = 0.0;
        com_mid[i][1] = 0.0;
        com_mid[i][2] = 0.0;
    }
    for(i=0; i<nbin; i++)
    {
        avgRDF[i] = 0;
    }
    //Open files
    FILE *traj, *isc_comtraj, *RDFfile;
    char *filename = malloc(sizeof(char)*35);
    char *filename3 = malloc(sizeof(char)*35);
    char *filename2 = malloc(sizeof(char)*35);
    char* str = malloc(sizeof(char)*35);

    sprintf(filename,"traj_Nbb%d_Nsc%d_f%d.xyz",Nbb, Nsc, f_branch);
    //sprintf(filename3,"isctraj_Nbb%d_Nsc%d_f%d.xyz",Nbb, Nsc, f_branch);
    //sprintf(filename2,"isc_comtraj_Nbb%d_Nsc%d_f%d.xyz",Nbb, Nsc, f_branch);
    sprintf(str, "RDF0_discrete.txt");

    traj = fopen(filename, "r"); 
    //isctraj = fopen(filename3, "w");
    //isc_comtraj = fopen(filename2, "w");
    RDFfile = fopen(str, "w");

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
        //fprintf(isctraj, "%d\n%d 1.500000\n", total_beadsinISC, time[t]);
        //fprintf(isc_comtraj, "%d\n%d 1.500000\n", Tot_ISC, time[t]);
   
    // Map ISC segment 
        // Iterate over chains
        for(c=0; c<nchain; c++)
        {
            bb_beg = c*(natom_perchain);
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
                    // comx += rx;
                    // comy += ry;
                    // comz += rz;
                    com[c*ISC_segment+i][0] += rx;
                    com[c*ISC_segment+i][1] += ry;
                    com[c*ISC_segment+i][2] += rz;
                    count++;
                    //fprintf(isctraj, "%d %lf %lf %lf\n", isc_type[c*ISC_segment+i], coords[bb_isc_beg+j][0], coords[bb_isc_beg+j][1], coords[bb_isc_beg+j][2]);       
                    
                    for(k=0; k<f_branch; k++)
                    {
                        sc_isc_beg = bb_beg+Nbb+((i*bb_perISC+j) * f_branch + k) * Nsc;
                        for(l=0; l<Nsc; l++)
                        {
                            rxo = coords[sc_isc_beg+l][0];
                            ryo = coords[sc_isc_beg+l][1];
                            rzo = coords[sc_isc_beg+l][2];
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
                            //fprintf(isctraj, "%d %lf %lf %lf\n", isc_type[c*ISC_segment+i], coords[sc_isc_beg+l][0], coords[sc_isc_beg+l][1], coords[sc_isc_beg+l][2]);
                            //printf("sc beg %d\n", sc_isc_beg);
                        }
                    }
                }
                //ISC coordinates as COM 
                com[c*ISC_segment+i][0] /= count;
                com[c*ISC_segment+i][1] /= count;
                com[c*ISC_segment+i][2] /= count;
                //fprintf(isc_comtraj, "%d 1 %lf %lf %lf\n", c, com[c*ISC_segment+i][0], com[c*ISC_segment+i][1], com[c*ISC_segment+i][2]);
                //fprintf(isc_comtraj, "%d %lf %lf %lf\n", c, comx, comy, comz);
            }                   
        }
        //Discretize the ISC segment
        for(c=0; c<nchain; c++)
        {
            for(i=0; i<(ISC_segment-1); i++)
            {
                com_mid[c*(ISC_segment-1)+i][0] = (com[c*ISC_segment+i][0] + com[c*ISC_segment+i+1][0])*0.5; 
                com_mid[c*(ISC_segment-1)+i][1] = (com[c*ISC_segment+i][1] + com[c*ISC_segment+i+1][1])*0.5;
                com_mid[c*(ISC_segment-1)+i][2] = (com[c*ISC_segment+i][2] + com[c*ISC_segment+i+1][2])*0.5;
            }
        }
        // // calculate the RDF
        for(i=0; i<nbin; i++)
        {
            histogram[i] = 0;
        }

        for (i = 0; i < Tot_ISC - 1; i++) 
        {
            for (j = i + 1; j < Tot_ISC; j++) 
            {
                rx1 = com[i][0];
                ry1 = com[i][1];
                rz1 = com[i][2];
                rx2 = com[j][0];
                ry2 = com[j][1];
                rz2 = com[j][2];
                dx = simbox_1dgetimage((rx2 - rx1), box_size); 
                dy = simbox_1dgetimage((ry2 - ry1), box_size);
                dz = simbox_1dgetimage((rz2 - rz1), box_size);
                r = sqrt(dx*dx + dy*dy + dz*dz);
                bin = (int)(r / dr);
                if (bin < nbin) 
                {
                    histogram[bin]++;
                }
            }
        }
        for (i = 0; i < (nchain*(ISC_segment-1)) - 1; i++) 
        {
            for (j = i + 1; j < (nchain*(ISC_segment-1)); j++) 
            {
                rx1 = com_mid[i][0];
                ry1 = com_mid[i][1];
                rz1 = com_mid[i][2];
                rx2 = com_mid[j][0];
                ry2 = com_mid[j][1];
                rz2 = com_mid[j][2];
                dx = simbox_1dgetimage((rx2 - rx1), box_size); 
                dy = simbox_1dgetimage((ry2 - ry1), box_size);
                dz = simbox_1dgetimage((rz2 - rz1), box_size);
                r = sqrt(dx*dx + dy*dy + dz*dz);
                bin = (int)(r / dr); 
                if (bin < nbin) 
                {
                    histogram[bin]++;
                    
                }
            }
        }
        for (i = 0; i < (nchain*(ISC_segment)); i++) 
        {
            for (j = 0; j < (nchain*(ISC_segment-1)); j++) 
            {
                if((int)i/ISC_segment ==(int)j/(ISC_segment-1)) continue; // && abs(i%(ISC_segment)-(int)j%(ISC_segment-1)) < 2) continue;
                else if(i!=j)
                {
                    rx1 = com[i][0];
                    ry1 = com[i][1];
                    rz1 = com[i][2];
                    rx2 = com_mid[j][0];
                    ry2 = com_mid[j][1];
                    rz2 = com_mid[j][2];
                    dx = simbox_1dgetimage((rx2 - rx1), box_size); 
                    dy = simbox_1dgetimage((ry2 - ry1), box_size);
                    dz = simbox_1dgetimage((rz2 - rz1), box_size);
                    r = sqrt(dx*dx + dy*dy + dz*dz);
                    bin = (int)(r / dr); 
                    if (bin < nbin) 
                    {
                        histogram[bin]++;
                        
                    }
                }
            }
        }

        for (i = 0; i < nbin; i++) 
        {
            tot_gr = nchain*(ISC_segment+ISC_segment-1);
            r_low = i * dr;
            r_high = (i + 1) * dr;
            shell_volume = 4*M_PI*r_low*r_low*dr;
            //shell_volume = (4.0 / 3.0) * M_PI * (r_high * r_high * r_high - r_low * r_low * r_low);
            number_density = tot_gr / (box_size * box_size * box_size);
            normalization = (double) tot_gr * shell_volume * number_density;
            gr = (double) 2*histogram[i] / normalization; //We didn't count overlap pairs. Multiply 2
            if(t > (tmax-50)) 
            {
                avgRDF[i] += gr;            
            }
        }
    }

    for (i=0; i<nbin; i++)
    {
        fprintf(RDFfile, "%lf %lf\n", i*dr, avgRDF[i]/49);
    }
    fclose(RDFfile);

    //return 0;

}
double simbox_1dgetimage(double rx, double bl)
{
    double res = 0;
    res = rx - bl * round(rx / bl);
    return res;
}


