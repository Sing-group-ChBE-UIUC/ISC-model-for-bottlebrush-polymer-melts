#include "header.h"
#include "math.h"
void calc_gr_ISC(double **coords, double *avgRDF_AA, double *avgRDF_BB, double *avgRDF_AB, int *particle_type, int Nbb1, int Nbb2, int Nsc1, int Nsc2, int f_branch1, int f_branch2, int nchain1, int nchain2, int natom_perchain1, int natom_perchain2, double box_size,int NISC1, int NISC2, int bbperISC1, int bbperISC2, int nbin)
{
    int i, j, k, l, c, count, valid_nbb;
    int natom_total1, natom_total2, Tot_ISC1, Tot_ISC2, Tot_midISC1, Tot_midISC2, excluded_bb1, excluded_bb2;
    int bb_beg, bb_isc_beg, sc_isc_beg, bin;
    double dx, dy, dz, dxo, dyo, dzo, rx, ry, rz, rxo, ryo, rzo, r;
    double rx1, ry1, rz1, rx2, ry2, rz2;
    double dr, shell_volume, number_density, normalization, gr, r_low, r_high;
    
    natom_total1 = nchain1 * natom_perchain1;
    natom_total2 = nchain2 * natom_perchain2;
    dr = box_size / (2 * nbin);
    Tot_ISC1 = nchain1 * NISC1;
    Tot_ISC2 = nchain2 * NISC2;
    Tot_midISC1 = nchain1 * (NISC1-1);
    Tot_midISC2 = nchain2 * (NISC2-1);
    excluded_bb1 = Nbb1 - bbperISC1*NISC1;
    excluded_bb2 = Nbb2 - bbperISC2*NISC2;
    
    double **com1 = (double **)calloc(Tot_ISC1, sizeof(double *));
    for(i=0; i<(Tot_ISC1); i++)
    {
        com1[i] = (double *)calloc(3, sizeof(double));
    }
    double **com_mid1 = (double **)calloc(Tot_midISC1, sizeof(double *));
    for(i=0; i<(Tot_midISC1); i++)
    {
        com_mid1[i] = (double *)calloc(3, sizeof(double));
    } 
    double **com2 = (double **)calloc(Tot_ISC2, sizeof(double *));
    for(i=0; i<(Tot_ISC2); i++)
    {
        com2[i] = (double *)calloc(3, sizeof(double));
    }
    double **com_mid2 = (double **)calloc(Tot_midISC2, sizeof(double *));
    for(i=0; i<(Tot_midISC2); i++)
    {
        com_mid2[i] = (double *)calloc(3, sizeof(double));
    } 
    
    int *histogram_AA = calloc(nbin, sizeof(int));
    int *histogram_BB = calloc(nbin, sizeof(int));
    int *histogram_AB = calloc(nbin, sizeof(int));

    //Initialize isc type and RDF
    for(i=0; i<(Tot_ISC1); i++)
    {
        com1[i][0] = 0.0;
        com1[i][1] = 0.0;
        com1[i][2] = 0.0;
    }
    for(i=0; i<(Tot_midISC1); i++)
    {
        com_mid1[i][0] = 0.0;
        com_mid1[i][1] = 0.0;
        com_mid1[i][2] = 0.0;
    }
    for(i=0; i<nbin; i++)
    {
        avgRDF_AA[i] = 0;
        avgRDF_BB[i] = 0;
        avgRDF_AB[i] = 0;
    }
    for(i=0; i<(Tot_ISC2); i++)
    {
        com2[i][0] = 0.0;
        com2[i][1] = 0.0;
        com2[i][2] = 0.0;
    }
    for(i=0; i<(Tot_midISC2); i++)
    {
        com_mid2[i][0] = 0.0;
        com_mid2[i][1] = 0.0;
        com_mid2[i][2] = 0.0;
    }
    
    // Iterate over chains nchain1
    for(c=0; c<(nchain1); c++)
    {
        bb_beg = c*(natom_perchain1);
        for(i=0; i<NISC1; i++)
        {
            // Backbone beads within ISC segment
            count=0;
            dx=0; dy=0; dz=0; dxo=0; dyo=0; dzo=0;
            bb_isc_beg = bb_beg + i*bbperISC1;

            // consider excluded backbone 
            if(excluded_bb1 <0 && i == (NISC1-1)) valid_nbb = bbperISC1 + excluded_bb1;
            else valid_nbb = bbperISC1;
            for(j=1; j<valid_nbb; j++)
            {
                rxo = coords[bb_isc_beg+j][0];
                ryo = coords[bb_isc_beg+j][1];
                rzo = coords[bb_isc_beg+j][2];
                dx = rxo-coords[bb_isc_beg][0];
                dy = ryo-coords[bb_isc_beg][1];
                dz = rzo-coords[bb_isc_beg][2];
                if(dx > box_size/2) rxo -= box_size; if(dx < -box_size/2) rxo += box_size;
                if(dy > box_size/2) ryo -= box_size; if(dy < -box_size/2) ryo += box_size;
                if(dz > box_size/2) rzo -= box_size; if(dz < -box_size/2) rzo += box_size;
                com1[c*NISC1+i][0] += rxo;
                com1[c*NISC1+i][1] += ryo;
                com1[c*NISC1+i][2] += rzo;
                count++;

                for(k=0; k<f_branch1; k++)
                {
                    sc_isc_beg = bb_beg+Nbb1+((i*bbperISC1+j) * f_branch1 + k) * Nsc1;
                    for(l=0; l<Nsc1; l++)
                    {
                        rxo = coords[sc_isc_beg+l][0];
                        ryo = coords[sc_isc_beg+l][1];
                        rzo = coords[sc_isc_beg+l][2];
                        dx = rxo-coords[bb_isc_beg][0];
                        dy = ryo-coords[bb_isc_beg][1];
                        dz = rzo-coords[bb_isc_beg][2];
                        if(dx > box_size/2) rxo -= box_size; if(dx < -box_size/2) rxo += box_size;
                        if(dy > box_size/2) ryo -= box_size; if(dy < -box_size/2) ryo += box_size;
                        if(dz > box_size/2) rzo -= box_size; if(dz < -box_size/2) rzo += box_size;
                        com1[c*NISC1+i][0] += rxo;
                        com1[c*NISC1+i][1] += ryo;
                        com1[c*NISC1+i][2] += rzo;
                        count++;
                    }
                }
            }
            //ISC coordinates as COM 
            com1[c*NISC1+i][0] /= count;
            com1[c*NISC1+i][0]= simbox_1dgetimage(com1[c*NISC1+i][0], box_size);
            com1[c*NISC1+i][1] /= count;
            com1[c*NISC1+i][1]= simbox_1dgetimage(com1[c*NISC1+i][1], box_size);
            com1[c*NISC1+i][2] /= count;
            com1[c*NISC1+i][2]= simbox_1dgetimage(com1[c*NISC1+i][2], box_size);
        }                   
    }
    // Iterate over chains nchain2
    for(c=0; c<(nchain2); c++)
    {
        bb_beg = nchain1*natom_perchain1+c*(natom_perchain2);
        // printf("bb_beg: %d\n", bb_beg);
        for(i=0; i<NISC2; i++)
        {
            // Backbone beads within ISC segment
            count=0;
            dx=0; dy=0; dz=0; dxo=0; dyo=0; dzo=0;
            bb_isc_beg = bb_beg + i*bbperISC2;

            // consider excluded backbone 
            if(excluded_bb2 <0 && i == (NISC2-1)) valid_nbb = bbperISC2 + excluded_bb2;
            else valid_nbb = bbperISC2;
            for(j=1; j<valid_nbb; j++)
            {
                rxo = coords[bb_isc_beg+j][0];
                ryo = coords[bb_isc_beg+j][1];
                rzo = coords[bb_isc_beg+j][2];
                dx = rxo-coords[bb_isc_beg][0];
                dy = ryo-coords[bb_isc_beg][1];
                dz = rzo-coords[bb_isc_beg][2];
                if(dx > box_size/2) rxo -= box_size; if(dx < -box_size/2) rxo += box_size;
                if(dy > box_size/2) ryo -= box_size; if(dy < -box_size/2) ryo += box_size;
                if(dz > box_size/2) rzo -= box_size; if(dz < -box_size/2) rzo += box_size;
                com2[c*NISC2+i][0] += rxo;
                com2[c*NISC2+i][1] += ryo;
                com2[c*NISC2+i][2] += rzo;
                count++;
                for(k=0; k<f_branch2; k++)
                {
                    sc_isc_beg = bb_beg+Nbb2+((i*bbperISC2+j) * f_branch2 + k) * Nsc2;
                    // printf("sc_isc_beg: %d\n", sc_isc_beg);
                    for(l=0; l<Nsc2; l++)
                    {
                        rxo = coords[sc_isc_beg+l][0];
                        ryo = coords[sc_isc_beg+l][1];
                        rzo = coords[sc_isc_beg+l][2];
                        dx = rxo-coords[bb_isc_beg][0];
                        dy = ryo-coords[bb_isc_beg][1];
                        dz = rzo-coords[bb_isc_beg][2];
                        if(dx > box_size/2) rxo -= box_size; if(dx < -box_size/2) rxo += box_size;
                        if(dy > box_size/2) ryo -= box_size; if(dy < -box_size/2) ryo += box_size;
                        if(dz > box_size/2) rzo -= box_size; if(dz < -box_size/2) rzo += box_size;
                        com2[c*NISC2+i][0] += rxo;
                        com2[c*NISC2+i][1] += ryo;
                        com2[c*NISC2+i][2] += rzo;
                        count++;
                    }
                }
            }
            //ISC coordinates as COM 
            com2[c*NISC2+i][0] /= count;
            com2[c*NISC2+i][0]= simbox_1dgetimage(com2[c*NISC2+i][0], box_size);
            com2[c*NISC2+i][1] /= count;
            com2[c*NISC2+i][1]= simbox_1dgetimage(com2[c*NISC2+i][1], box_size);
            com2[c*NISC2+i][2] /= count;
            com2[c*NISC2+i][2]= simbox_1dgetimage(com2[c*NISC2+i][2], box_size);
        }                   
    }

    // Distretize mid ISC beads nchain1
    for(c=0; c<nchain1; c++)
    {
        for(i=0; i<(NISC1-1); i++)
        {
            dx = com1[c*NISC1+i][0]-com1[c*NISC1+i+1][0];
            dy = com1[c*NISC1+i][1]-com1[c*NISC1+i+1][1];
            dz = com1[c*NISC1+i][2]-com1[c*NISC1+i+1][2];
            com_mid1[c*(NISC1-1)+i][0] = (com1[c*NISC1+i][0] + com1[c*NISC1+i+1][0])*0.5;
            if(dx > box_size/2) com_mid1[c*(NISC1-1)+i][0] = 0.5*((com1[c*NISC1+i][0] - box_size)+com1[c*NISC1+i+1][0]);
            if(dx < -box_size/2) com_mid1[c*(NISC1-1)+i][0] = 0.5*((com1[c*NISC1+i][0] + box_size)+com1[c*NISC1+i+1][0]);
            com_mid1[c*(NISC1-1)+i][1] = (com1[c*NISC1+i][1] + com1[c*NISC1+i+1][1])*0.5;
            if(dy > box_size/2) com_mid1[c*(NISC1-1)+i][1] = 0.5*((com1[c*NISC1+i][1] - box_size)+com1[c*NISC1+i+1][1]);
            if(dy < -box_size/2) com_mid1[c*(NISC1-1)+i][1] = 0.5*((com1[c*NISC1+i][1] + box_size)+com1[c*NISC1+i+1][1]);
            com_mid1[c*(NISC1-1)+i][2] = (com1[c*NISC1+i][2] + com1[c*NISC1+i+1][2])*0.5;
            if(dz > box_size/2) com_mid1[c*(NISC1-1)+i][2] = 0.5*((com1[c*NISC1+i][2] - box_size)+com1[c*NISC1+i+1][2]);
            if(dz < -box_size/2) com_mid1[c*(NISC1-1)+i][2] = 0.5*((com1[c*NISC1+i][2] + box_size)+com1[c*NISC1+i+1][2]);
        }
    }
    // Distretize mid ISC beads nchain2
    for(c=0; c<nchain2; c++)
    {
        for(i=0; i<(NISC2-1); i++)
        {
            dx = com2[c*NISC2+i][0]-com2[c*NISC2+i+1][0];
            dy = com2[c*NISC2+i][1]-com2[c*NISC2+i+1][1];
            dz = com2[c*NISC2+i][2]-com2[c*NISC2+i+1][2];
            com_mid2[c*(NISC2-1)+i][0] = (com2[c*NISC2+i][0] + com2[c*NISC2+i+1][0])*0.5;
            if(dx > box_size/2) com_mid2[c*(NISC2-1)+i][0] = 0.5*((com2[c*NISC2+i][0] - box_size)+com2[c*NISC2+i+1][0]);
            if(dx < -box_size/2) com_mid2[c*(NISC2-1)+i][0] = 0.5*((com2[c*NISC2+i][0] + box_size)+com2[c*NISC2+i+1][0]);
            com_mid2[c*(NISC2-1)+i][1] = (com2[c*NISC2+i][1] + com2[c*NISC2+i+1][1])*0.5;
            if(dy > box_size/2) com_mid2[c*(NISC2-1)+i][1] = 0.5*((com2[c*NISC2+i][1] - box_size)+com2[c*NISC2+i+1][1]);
            if(dy < -box_size/2) com_mid2[c*(NISC2-1)+i][1] = 0.5*((com2[c*NISC2+i][1] + box_size)+com2[c*NISC2+i+1][1]);
            com_mid2[c*(NISC2-1)+i][2] = (com2[c*NISC2+i][2] + com2[c*NISC2+i+1][2])*0.5;
            if(dz > box_size/2) com_mid2[c*(NISC2-1)+i][2] = 0.5*((com2[c*NISC2+i][2] - box_size)+com2[c*NISC2+i+1][2]);
            if(dz < -box_size/2) com_mid2[c*(NISC2-1)+i][2] = 0.5*((com2[c*NISC2+i][2] + box_size)+com2[c*NISC2+i+1][2]);
        }
    }
    
    // Calculate the RDF
    for(i=0; i<nbin; i++)
    {
        histogram_AA[i] = 0;
        histogram_BB[i] = 0;
        histogram_AB[i] = 0;
    }
 

    // RDF calculation
    // // A com-com
    // gr_prep(0, NISC1, NISC1, Tot_ISC1, Tot_ISC1, com1, com1, box_size, nbin, histogram_AA);
    // // A mid-mid
    // gr_prep(0, NISC1-1, NISC1-1, Tot_midISC1, Tot_midISC1, com_mid1, com_mid1, box_size, nbin, histogram_AA);
    // // A mid-com
    // gr_prep(2, NISC1, NISC1-1, Tot_ISC1, Tot_midISC1, com1, com_mid1, box_size, nbin, histogram_AA);
    // // B com-com
    // gr_prep(0, NISC2, NISC2, Tot_ISC2, Tot_ISC2, com2, com2, box_size, nbin, histogram_BB);
    // // B mid-mid
    // gr_prep(0, NISC2-1, NISC2-1, Tot_midISC2, Tot_midISC2, com_mid2, com_mid2, box_size, nbin, histogram_BB);
    // // B mid-com
    // gr_prep(2, NISC2, NISC2-1, Tot_ISC2, Tot_midISC2, com2, com_mid2, box_size, nbin, histogram_BB);
    // // A com-B com
    // gr_prep(1, NISC1, NISC2, Tot_ISC1, Tot_ISC2, com1, com2, box_size, nbin, histogram_AB);
    // // A com-B mid
    // gr_prep(1, NISC1, NISC2-1, Tot_ISC1, Tot_midISC2, com1, com_mid2, box_size, nbin, histogram_AB);
    // // A mid-B com
    // gr_prep(1, NISC1-1, NISC2, Tot_midISC1, Tot_ISC2, com_mid1, com2, box_size, nbin, histogram_AB);
    // // A mid-B mid
    // gr_prep(1, NISC1-1, NISC2-1, Tot_midISC1, Tot_midISC2, com_mid1, com_mid2, box_size, nbin, histogram_AB);

    // RDF original indexes
    for (i = 0; i < Tot_ISC1 - 1; i++) 
    {
        for (j = i + 1; j < Tot_ISC1; j++) 
        {
            rx1 = com1[i][0];
            ry1 = com1[i][1];
            rz1 = com1[i][2];
            rx2 = com1[j][0];
            ry2 = com1[j][1];
            rz2 = com1[j][2];
            dx = rx2-rx1; dy = ry2-ry1; dz = rz2-rz1;
            if(dx > box_size/2) dx -= box_size; if(dx < -box_size/2) dx += box_size;
            if(dy > box_size/2) dy -= box_size; if(dy < -box_size/2) dy += box_size;
            if(dz > box_size/2) dz -= box_size; if(dz < -box_size/2) dz += box_size;
            r = sqrt(dx*dx + dy*dy + dz*dz);
            bin = (int)(r / dr);
            if (bin < nbin) histogram_AA[bin]++;
        }
    }
    // RDF mid ISC indexes
    for (i = 0; i < Tot_midISC1-1; i++) 
    {
        for (j = i+1; j < Tot_midISC1; j++) 
        {
            rx1 = com_mid1[i][0];
            ry1 = com_mid1[i][1];
            rz1 = com_mid1[i][2];
            rx2 = com_mid1[j][0];
            ry2 = com_mid1[j][1];
            rz2 = com_mid1[j][2];
            dx = rx2-rx1; dy = ry2-ry1; dz = rz2-rz1;
            if(dx > box_size/2) dx -= box_size; if(dx < -box_size/2) dx += box_size;
            if(dy > box_size/2) dy -= box_size; if(dy < -box_size/2) dy += box_size;
            if(dz > box_size/2) dz -= box_size; if(dz < -box_size/2) dz += box_size;
            r = sqrt(dx*dx + dy*dy + dz*dz);
            bin = (int)(r / dr);
            if (bin < nbin) histogram_AA[bin]++;
        }
    }
    // RDF between original and mid ISC
    for (i = 0; i < Tot_ISC1; i++) 
    {
        for (j = 0; j < Tot_midISC1; j++) 
        {
            if((int)i/(NISC1) == (int)j/(NISC1-1) && (abs((int)i%(NISC1)-(int)j%(NISC1-1)) < 2)) continue;
            else
            {
                rx1 = com1[i][0];
                ry1 = com1[i][1];
                rz1 = com1[i][2];
                rx2 = com_mid1[j][0];
                ry2 = com_mid1[j][1];
                rz2 = com_mid1[j][2];
                dx = rx2-rx1; dy = ry2-ry1; dz = rz2-rz1;   
                if(dx > box_size/2) dx -= box_size; if(dx < -box_size/2) dx += box_size;
                if(dy > box_size/2) dy -= box_size; if(dy < -box_size/2) dy += box_size;
                if(dz > box_size/2) dz -= box_size; if(dz < -box_size/2) dz += box_size;
                r = sqrt(dx*dx + dy*dy + dz*dz);
                bin = (int)(r / dr);
                if (bin < nbin) histogram_AA[bin]++;
            }
        }
    }
    for (i = 0; i < Tot_ISC2 - 1; i++) 
    {
        for (j = i + 1; j < Tot_ISC2; j++) 
        {
            rx1 = com2[i][0];
            ry1 = com2[i][1];
            rz1 = com2[i][2];
            rx2 = com2[j][0];
            ry2 = com2[j][1];
            rz2 = com2[j][2];
            dx = rx2-rx1; dy = ry2-ry1; dz = rz2-rz1;
            if(dx > box_size/2) dx -= box_size; if(dx < -box_size/2) dx += box_size;
            if(dy > box_size/2) dy -= box_size; if(dy < -box_size/2) dy += box_size;
            if(dz > box_size/2) dz -= box_size; if(dz < -box_size/2) dz += box_size;
            r = sqrt(dx*dx + dy*dy + dz*dz);
            bin = (int)(r / dr);
            if (bin < nbin) histogram_BB[bin]++;
        }
    }
    // RDF mid ISC indexes
    for (i = 0; i < Tot_midISC2-1; i++) 
    {
        for (j = i+1; j < Tot_midISC2; j++) 
        {
            rx1 = com_mid2[i][0];
            ry1 = com_mid2[i][1];
            rz1 = com_mid2[i][2];
            rx2 = com_mid2[j][0];
            ry2 = com_mid2[j][1];
            rz2 = com_mid2[j][2];
            dx = rx2-rx1; dy = ry2-ry1; dz = rz2-rz1;
            if(dx > box_size/2) dx -= box_size; if(dx < -box_size/2) dx += box_size;
            if(dy > box_size/2) dy -= box_size; if(dy < -box_size/2) dy += box_size;
            if(dz > box_size/2) dz -= box_size; if(dz < -box_size/2) dz += box_size;
            r = sqrt(dx*dx + dy*dy + dz*dz);
            bin = (int)(r / dr);
            if (bin < nbin) histogram_BB[bin]++;
        }
    }
    // RDF between original and mid ISC
    for (i = 0; i < Tot_ISC2; i++) 
    {
        for (j = 0; j < Tot_midISC2; j++) 
        {
            if((int)i/(NISC2) == (int)j/(NISC2-1) && (abs((int)i%(NISC2)-(int)j%(NISC2-1)) < 2)) continue;
            else
            {
                rx1 = com2[i][0];
                ry1 = com2[i][1];
                rz1 = com2[i][2];
                rx2 = com_mid2[j][0];
                ry2 = com_mid2[j][1];
                rz2 = com_mid2[j][2];
                dx = rx2-rx1; dy = ry2-ry1; dz = rz2-rz1;   
                if(dx > box_size/2) dx -= box_size; if(dx < -box_size/2) dx += box_size;
                if(dy > box_size/2) dy -= box_size; if(dy < -box_size/2) dy += box_size;
                if(dz > box_size/2) dz -= box_size; if(dz < -box_size/2) dz += box_size;
                r = sqrt(dx*dx + dy*dy + dz*dz);
                bin = (int)(r / dr);
                if (bin < nbin) histogram_BB[bin]++;
            }
        }
    }

    for (i = 0; i < Tot_ISC1; i++) 
    {
        for (j = 0; j < Tot_ISC2; j++) 
        {
            rx1 = com1[i][0];
            ry1 = com1[i][1];
            rz1 = com1[i][2];
            rx2 = com2[j][0];
            ry2 = com2[j][1];
            rz2 = com2[j][2];
            dx = rx2-rx1; dy = ry2-ry1; dz = rz2-rz1;
            if(dx > box_size/2) dx -= box_size; if(dx < -box_size/2) dx += box_size;
            if(dy > box_size/2) dy -= box_size; if(dy < -box_size/2) dy += box_size;
            if(dz > box_size/2) dz -= box_size; if(dz < -box_size/2) dz += box_size;
            r = sqrt(dx*dx + dy*dy + dz*dz);
            bin = (int)(r / dr);
            if (bin < nbin) histogram_AB[bin]++;
        }
    }
    for (i = 0; i < Tot_ISC1; i++) 
    {
        for (j = 0; j < Tot_midISC2; j++) 
        {
            rx1 = com1[i][0];
            ry1 = com1[i][1];
            rz1 = com1[i][2];
            rx2 = com_mid2[j][0];
            ry2 = com_mid2[j][1];
            rz2 = com_mid2[j][2];
            dx = rx2-rx1; dy = ry2-ry1; dz = rz2-rz1;
            if(dx > box_size/2) dx -= box_size; if(dx < -box_size/2) dx += box_size;
            if(dy > box_size/2) dy -= box_size; if(dy < -box_size/2) dy += box_size;
            if(dz > box_size/2) dz -= box_size; if(dz < -box_size/2) dz += box_size;
            r = sqrt(dx*dx + dy*dy + dz*dz);
            bin = (int)(r / dr);
            if (bin < nbin) histogram_AB[bin]++;
        }
    }
    for (i = 0; i < Tot_midISC1; i++) 
    {
        for (j = 0; j < Tot_midISC2; j++) 
        {
            rx1 = com_mid1[i][0];
            ry1 = com_mid1[i][1];
            rz1 = com_mid1[i][2];
            rx2 = com_mid2[j][0];
            ry2 = com_mid2[j][1];
            rz2 = com_mid2[j][2];
            dx = rx2-rx1; dy = ry2-ry1; dz = rz2-rz1;
            if(dx > box_size/2) dx -= box_size; if(dx < -box_size/2) dx += box_size;
            if(dy > box_size/2) dy -= box_size; if(dy < -box_size/2) dy += box_size;
            if(dz > box_size/2) dz -= box_size; if(dz < -box_size/2) dz += box_size;
            r = sqrt(dx*dx + dy*dy + dz*dz);
            bin = (int)(r / dr);
            if (bin < nbin) histogram_AB[bin]++;
        }
    }
    for (i = 0; i < Tot_midISC1; i++) 
    {
        for (j = 0; j < Tot_ISC2; j++) 
        {
            rx1 = com_mid1[i][0];
            ry1 = com_mid1[i][1];
            rz1 = com_mid1[i][2];
            rx2 = com2[j][0];
            ry2 = com2[j][1];
            rz2 = com2[j][2];
            dx = rx2-rx1; dy = ry2-ry1; dz = rz2-rz1;
            if(dx > box_size/2) dx -= box_size; if(dx < -box_size/2) dx += box_size;
            if(dy > box_size/2) dy -= box_size; if(dy < -box_size/2) dy += box_size;
            if(dz > box_size/2) dz -= box_size; if(dz < -box_size/2) dz += box_size;
            r = sqrt(dx*dx + dy*dy + dz*dz);
            bin = (int)(r / dr);
            if (bin < nbin) histogram_AB[bin]++;
        }
    }

    // AA Normalize RDF
    for (i = 0; i < nbin; i++) 
    {
        r_low = i * dr;
        r_high = (i + 1) * dr;
        shell_volume = 4*M_PI*r_low*r_low*dr;
        number_density = (Tot_ISC1+Tot_midISC1) / (box_size * box_size * box_size);
        normalization = (double) (Tot_ISC1+Tot_midISC1) * shell_volume * number_density;
        gr = (double) 2*histogram_AA[i]/ normalization ; //We didn't count overlap pairs. Multiply 2
        avgRDF_AA[i] = gr;       
    }

    // BB Normalize RDF
    for (i = 0; i < nbin; i++) 
    {
        r_low = i * dr;
        r_high = (i + 1) * dr;
        shell_volume = 4*M_PI*r_low*r_low*dr;
        number_density = (Tot_ISC2+Tot_midISC2) / (box_size * box_size * box_size);
        normalization = (double) (Tot_ISC2+Tot_midISC2) * shell_volume * number_density;
        gr = (double) 2*histogram_BB[i] / normalization; //We didn't count overlap pairs. Multiply 2
        avgRDF_BB[i] = gr;            
    }

    //AB Normalize RDF
    for (i = 0; i < nbin; i++) 
    {
        r_low = i * dr;
        r_high = (i + 1) * dr;
        shell_volume = 4*M_PI*r_low*r_low*dr;
        number_density = (Tot_ISC2+Tot_midISC2) / (box_size * box_size * box_size);
        normalization = (double) (Tot_ISC1+Tot_midISC1) * shell_volume * number_density;
        gr = (double) histogram_AB[i] / normalization;
        avgRDF_AB[i] = gr;           
    }

    // Free memory 
    for(i=0; i<(Tot_ISC1); i++)
    {
        free(com1[i]);
    }
    free(com1);
    for(i=0; i<(Tot_midISC1); i++)
    {
        free(com_mid1[i]);
    }
    free(com_mid1);
    for(i=0; i<(Tot_ISC2); i++)
    {
        free(com2[i]);
    }
    free(com2);
    for(i=0; i<(Tot_midISC2); i++)
    {
        free(com_mid2[i]);
    }
    free(com_mid2);
    free(histogram_AA);
    free(histogram_BB);
    free(histogram_AB);
}

