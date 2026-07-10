/////////////////////////////
// Created by Haisu Kang 05/19/2024
// Cautious that some of the lines are assuming 
//    nchain1 = nchain2, natom_perchain1 = natom_perchain2, Nbb1 = Nbb2, Nsc1 = Nsc2, f_branch1 = f_branch2
////////////////////////////////

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

double simbox_1dgetimage(double rx, double bl);

int main (int argc, const char * argv[])
{
    //Variables
    int i, j, k, l, c, t, nbin, bin, bb_beg, bb_isc_beg, sc_isc_beg, particle_type_index_i, particle_type_index_j, valid_nbb, excluded_bb, count, Nbb, Nsc, f_branch;
    int nchain1, nchain2, natom_perchain1, natom_perchain2, natom_total, nbonds_total, nint, tmax, nsites, maxsite_1d, ISC_segment, Tot_ISC1, Tot_ISC2, bb_perISC, total_bbinISC;
    double density, rho_norm, bondlen, chi, kappa, kspring, box_size, grid_size, box_volume, dr, dx, dy, dz, dxo, dyo, dzo, rx, ry, rz, rxo, ryo, rzo, comx, comy, comz, r, r_low, r_high, shell_volume, number_density, normalization, gr;
    double rx1, ry1, rz1, rx2, ry2, rz2;
    //Read output file
    FILE *outputfile;
    outputfile = fopen("output.txt", "r");
    fscanf(outputfile, "nchain1 = %d\n", &nchain1);
    fscanf(outputfile, "nchain2 = %d\n", &nchain2);
    fscanf(outputfile, "natom_perchain1 = %d\n", &natom_perchain1);
    fscanf(outputfile, "natom_perchain2 = %d\n", &natom_perchain2);
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
    box_size = maxsite_1d * grid_size;
    box_volume = box_size * box_size * box_size;
    tmax = 20;
    nbin = 500;
    dr = box_size / (2.0 * nbin);
    Nsc = 1;
    Nbb = 11;
    f_branch = 1;

    //ISC parameters
    ISC_segment =  atoi(argv[1]);
    Tot_ISC1 = ISC_segment*nchain1;
    Tot_ISC2 = ISC_segment*nchain2;
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
    double **com = (double **)calloc(Tot_ISC1+Tot_ISC2, sizeof(double *));
    for(i=0; i<(Tot_ISC1+Tot_ISC2); i++)
    {
        com[i] = (double *)calloc(3, sizeof(double));
    }
    double **com_mid = (double **)calloc(2*nchain1*(ISC_segment-1), sizeof(double *));
    for(i=0; i<(2*nchain1*(ISC_segment-1)); i++)
    {
        com_mid[i] = (double *)calloc(3, sizeof(double));
    }
    int *particle_type = calloc(natom_total, sizeof(int));
    int *time = calloc(tmax, sizeof(int));
    int *histogram_AA = calloc(nbin, sizeof(int));
    int *histogram_AB = calloc(nbin, sizeof(int));
    double *avgRDF_AA = calloc(nbin, sizeof(double));
    double *avgRDF_AB = calloc(nbin, sizeof(double));

    //Initialize isc type and RDF
    for(i=0; i<(Tot_ISC1+Tot_ISC2); i++)
    {
        com[i][0] = 0.0;
        com[i][1] = 0.0;
        com[i][2] = 0.0;
    }
    for (i=0; i<(2*nchain1*(ISC_segment-1)); i++)
    {
        com_mid[i][0] = 0.0;
        com_mid[i][1] = 0.0;
        com_mid[i][2] = 0.0;
    }
    for(i=0; i<nbin; i++)
    {
        avgRDF_AA[i] = 0;
        avgRDF_AB[i] = 0;
    }
    //Open files
    FILE *traj, *RDFfile_AA, *RDFfile_AB;
    char *filename = malloc(sizeof(char)*35);
    char* str = malloc(sizeof(char)*35);
    char* str1 = malloc(sizeof(char)*35);

    sprintf(filename,"traj_blend.xyz");
    sprintf(str, "RDF_AA.txt");
    sprintf(str1, "RDF_AB.txt");

    traj = fopen(filename, "r"); 
    RDFfile_AA = fopen(str, "w");
    RDFfile_AB = fopen(str1, "w");


    // Read the trajectory file
    for (t=0; t<tmax; t++)
    {
        // read the coordinates of the particles
        fscanf(traj, "%d\n%d 1.500000\n", &natom_total, &time[t]);
        for (i = 0; i < natom_total; ++i)
        {
            fscanf(traj, "%d %lf %lf %lf\n", &particle_type[i], &coords[i][0], &coords[i][1], &coords[i][2]);
        }
   
    // Map ISC segment 
        // Iterate over chains
        for(c=0; c<(nchain1+nchain2); c++)
        {
            bb_beg = c*(natom_perchain1);
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
                        }
                    }
                }
                //ISC coordinates as COM 
                com[c*ISC_segment+i][0] /= count;
                com[c*ISC_segment+i][1] /= count;
                com[c*ISC_segment+i][2] /= count;
            }                   
        }
        //Discretize the ISC segment A
        for(c=0; c<nchain1; c++)
        {
            for(i=0; i<(ISC_segment-1); i++)
            {
                com_mid[c*(ISC_segment-1)+i][0] = (com[c*ISC_segment+i][0] + com[c*ISC_segment+i+1][0])*0.5; 
                com_mid[c*(ISC_segment-1)+i][1] = (com[c*ISC_segment+i][1] + com[c*ISC_segment+i+1][1])*0.5;
                com_mid[c*(ISC_segment-1)+i][2] = (com[c*ISC_segment+i][2] + com[c*ISC_segment+i+1][2])*0.5;
            }
        }
        //Distretize the ISC segment B
        for(c=nchain1; c<(nchain1+nchain2); c++)
        {
            for(i=0; i<(ISC_segment-1); i++)
            {
                com_mid[c*(ISC_segment-1)+i][0] = (com[c*ISC_segment+i][0] + com[c*ISC_segment+i+1][0])*0.5; 
                com_mid[c*(ISC_segment-1)+i][1] = (com[c*ISC_segment+i][1] + com[c*ISC_segment+i+1][1])*0.5;
                com_mid[c*(ISC_segment-1)+i][2] = (com[c*ISC_segment+i][2] + com[c*ISC_segment+i+1][2])*0.5;
            }
        }
        
    // Calculate the RDF
        for(i=0; i<nbin; i++)
        {
            histogram_AA[i] = 0;
            histogram_AB[i] = 0;
        }
        // AA RDF original index
        for (i = 0; i < Tot_ISC1 - 1; i++) 
        {
            for (j = i + 1; j < Tot_ISC1; j++) 
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
                    histogram_AA[bin]++;
                }
                //printf("histogram_AA[%d] = %d\n", bin, histogram_AA[bin]);
            }
        }
        // AA RDF discretized index (midpoints)
        for (i=0; i<(nchain1*(ISC_segment-1))-1; i++)
        {
            for(j=i+1; j<(nchain1*(ISC_segment-1)); j++)
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
                    histogram_AA[bin]++;
                }
                //printf("histogram_BB[%d] = %d\n", bin, histogram_BB[bin]);
            }
        }
        // AA RDF between original and discretized index
        for (i = 0; i < (nchain1*(ISC_segment)); i++) 
        {
            for (j = 0; j < (nchain1*(ISC_segment-1)); j++) 
            {
                if((int)i/ISC_segment ==(int)j/(ISC_segment-1)) continue; // && abs(i%(ISC_segment)-(int)j%(ISC_segment-1)) < 2) continue;
                else 
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
                        histogram_AA[bin]++;
                        
                    }
                }
            }
        }

        // AB RDF original index
        for (i = 0; i < Tot_ISC1; i++) 
        {
            for (j = 0; j < Tot_ISC2; j++) 
            {
                particle_type_index_i=(int)(i/ISC_segment)*natom_perchain1+(int)(i%ISC_segment)*bb_perISC+1;
                particle_type_index_j=(int)((j+Tot_ISC1)/ISC_segment)*natom_perchain1+(int)((j+Tot_ISC1)%ISC_segment)*bb_perISC+1;
                //printf("particle_type_index_i = %d, particle_type_index_j = %d\n", particle_type[particle_type_index_i], particle_type[particle_type_index_j]);
                rx1 = com[i][0];
                ry1 = com[i][1];
                rz1 = com[i][2];
                rx2 = com[Tot_ISC1+j][0];
                ry2 = com[Tot_ISC1+j][1];
                rz2 = com[Tot_ISC1+j][2];
                dx = simbox_1dgetimage((rx2 - rx1), box_size); 
                dy = simbox_1dgetimage((ry2 - ry1), box_size);
                dz = simbox_1dgetimage((rz2 - rz1), box_size);
                r = sqrt(dx*dx + dy*dy + dz*dz);
                bin = (int)(r / dr);
                if (bin < nbin) 
                {
                    histogram_AB[bin]++;
                }

            }
        }
        // AB RDF discretized index (midpoints)
        for (i = 0; i < (nchain1*(ISC_segment-1)); i++) 
        {
            for (j = 0; j < (nchain1*(ISC_segment-1)); j++) 
            {
                particle_type_index_i=(int)(i/(ISC_segment-1))*natom_perchain1+(int)(i%(ISC_segment-1))*bb_perISC+1;
                particle_type_index_j=(int)((nchain1*(ISC_segment-1)+j)/(ISC_segment-1))*natom_perchain1+(int)((nchain1*(ISC_segment-1)+j)%(ISC_segment-1))*bb_perISC+1;
                //printf("particle_type_index_i = %d, particle_type_index_j = %d\n", particle_type[particle_type_index_i], particle_type[particle_type_index_j]);
                rx1 = com_mid[i][0];
                ry1 = com_mid[i][1];
                rz1 = com_mid[i][2];
                rx2 = com_mid[nchain1*(ISC_segment-1)+j][0];
                ry2 = com_mid[nchain1*(ISC_segment-1)+j][1];
                rz2 = com_mid[nchain1*(ISC_segment-1)+j][2];
                dx = simbox_1dgetimage((rx2 - rx1), box_size); 
                dy = simbox_1dgetimage((ry2 - ry1), box_size);
                dz = simbox_1dgetimage((rz2 - rz1), box_size);
                r = sqrt(dx*dx + dy*dy + dz*dz);
                bin = (int)(r / dr);
                if (bin < nbin) 
                {
                    histogram_AB[bin]++;
                }

            }
        }
        // AB RDF between original and discretized index
        for (i = 0; i < (Tot_ISC1+Tot_ISC2); i++) 
        {
            for (j = 0; j < (2*nchain1*(ISC_segment-1)); j++) 
            {
                particle_type_index_i=(int)(i/ISC_segment)*natom_perchain1+(int)(i%ISC_segment)*bb_perISC+1;
                particle_type_index_j=(int)(j/(ISC_segment-1))*natom_perchain1+(int)(j%(ISC_segment-1))*bb_perISC+1;
                if(particle_type[particle_type_index_i] != particle_type[particle_type_index_j])
                {
                    //printf("particle_type_index_i = %d, particle_type_index_j = %d\n", particle_type[particle_type_index_i], particle_type[particle_type_index_j]);
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
                        histogram_AB[bin]++;
                    }
                }
            }
        }

        for (i = 0; i < nbin; i++) 
        {
            r_low = i * dr;
            r_high = (i + 1) * dr;
            shell_volume = 4*M_PI*r_low*r_low*dr;
            number_density = (nchain1*(2*ISC_segment-1)) / (box_size * box_size * box_size);
            normalization = (double) (nchain1*(2*ISC_segment-1)) * shell_volume * number_density;
            gr = (double) 2*histogram_AA[i]/ normalization ; //We didn't count overlap pairs. Multiply 2
            if(t > 1) 
            {
                avgRDF_AA[i] += gr;            
            }
        }

        for (i = 0; i < nbin; i++) 
        {
            r_low = i * dr;
            r_high = (i + 1) * dr;
            shell_volume = 4*M_PI*r_low*r_low*dr;
            number_density = (nchain1*(2*ISC_segment-1)) / (box_size * box_size * box_size);
            normalization = (double) ((nchain1*(2*ISC_segment-1)))*shell_volume * number_density;
            gr = (double) histogram_AB[i] / normalization; //We didn't count overlap pairs. Multiply 2
            if(t > 1) 
            {
                avgRDF_AB[i] += gr;            
            }
        }
    }

    for (i=0; i<nbin; i++)
    {
        fprintf(RDFfile_AA, "%lf %lf\n", i*dr, avgRDF_AA[i]/19);
    }
    fclose(RDFfile_AA);
    for (i=0; i<nbin; i++)
    {
        fprintf(RDFfile_AB, "%lf %lf\n", i*dr, avgRDF_AB[i]/19);
    }
    fclose(RDFfile_AB);
}



double simbox_1dgetimage(double rx, double bl)
{
    double res = 0;
    res = rx - bl * round(rx / bl);
    return res;
}


