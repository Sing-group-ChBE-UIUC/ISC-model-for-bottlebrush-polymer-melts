#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

int main (int argc, const char * argv[])
{
    int tot_idx, idx, x_idx, y_idx, z_idx, i, j, k, t, q, natom, tmax, tcount, q_index, dum, id, nq, num_bins;
    double  lx, ly, lz, x1, x2, y1, y2, z1, z2, rx, ry, rz, qmin, qmax, qbin, qxi, qyi, qzi, qx, qy, qz, del, oned_del;
    double real_sum, imag_sum, r, s_q_value, qmag;
    
    // Parameters
    qbin = 5;
    qmax = 0.6;
    qmin = 0.0001;
    num_bins = 100;
    nq = qbin*10*qbin*10*qbin*2;
    del = (double)((qmax-qmin)/num_bins);
    oned_del = (double)((qmax-qmin)/nq);
    tmax = 100;
    printf("del = %lf\n", del);
    printf("nq = %d\n", nq);

    FILE *traj = fopen(argv[1], "r");
    fscanf(traj,"ITEM: TIMESTEP\n%d\n", &dum);
    fscanf(traj,"ITEM: NUMBER OF ATOMS\n%d\n", &natom);
    // fscanf(traj,"ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
    fscanf(traj,"ITEM: BOX BOUNDS pp pp pp\n");
    fscanf(traj,"%lf %lf\n", &x1, &lx);
    fscanf(traj,"%lf %lf\n", &y1, &ly);
    fscanf(traj,"%lf %lf\n", &z1, &lz);
    // fscanf(traj,"%lf %lf %lf\n", &x1, &lx, &x2);
    // fscanf(traj,"%lf %lf %lf\n", &y1, &ly, &y2);
    // fscanf(traj,"%lf %lf %lf\n", &z1, &lz, &z2);
    // fscanf(traj,"ITEM: ATOMS id mol type x y z\n");
    fscanf(traj,"ITEM: ATOMS idtype x y z\n");
    fclose(traj);

    // Memory allocation
    int *timestep = calloc(tmax, sizeof(int));
    int *mol = calloc(natom, sizeof(int));
    int *type = calloc(natom, sizeof(int));
    double **x;
    x = (double**)calloc(tmax, sizeof(double *));
    for(i = 0; i < tmax; i++) x[i] = (double*)calloc(natom, sizeof(double));
    double **y;
    y = (double**)calloc(tmax, sizeof(double *));
    for(i = 0; i < tmax; i++) y[i] = (double*)calloc(natom, sizeof(double));
    double **z;
    z = (double**)calloc(tmax, sizeof(double *));
    for(i = 0; i < tmax; i++) z[i] = (double*)calloc(natom, sizeof(double));
    double *s_q = calloc(num_bins, sizeof(double));
    double *q_list = calloc(num_bins, sizeof(double));
    double *q_real = calloc(num_bins, sizeof(double));
    double *qdotr = calloc(lx*ly*lz, sizeof(double));
    double *xyz_grid_center = calloc(lx*ly*lz*3, sizeof(double));
    int *count_grid = calloc(lx*ly*lz, sizeof(int));
    int *count_bin = calloc(lx*ly*lz, sizeof(int));
    double *dens_grid = calloc(lx*ly*lz, sizeof(double));

    
    // Read dump file
    traj = fopen(argv[1], "r");
    for(t=0; t<tmax; t++)
    {
        fscanf(traj,"ITEM: TIMESTEP\n%d\n", &dum);
        fscanf(traj,"ITEM: NUMBER OF ATOMS\n%d\n", &natom);
        // fscanf(traj,"ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
        fscanf(traj,"ITEM: BOX BOUNDS pp pp pp\n");
        fscanf(traj,"%lf %lf\n", &x1, &lx);
        fscanf(traj,"%lf %lf\n", &y1, &ly);
        fscanf(traj,"%lf %lf\n", &z1, &lz);
        // fscanf(traj,"%lf %lf %lf\n", &x1, &lx, &x2);
        // fscanf(traj,"%lf %lf %lf\n", &y1, &ly, &y2);
        // fscanf(traj,"%lf %lf %lf\n", &z1, &lz, &z2);
        // fscanf(traj,"ITEM: ATOMS id mol type x y z\n");
        fscanf(traj,"ITEM: ATOMS idtype x y z\n");
        for (i = 0; i < natom; ++i) {
            fscanf(traj, "%d %d %lf %lf %lf\n", &id, &type[i], &x[t][i], &y[t][i], &z[t][i]);
        }

    }
    fclose(traj);
    lx=100; ly=100; lz=20;
   // Initialize s(q)
    for(q=0; q<num_bins; q++)
    {
        s_q[q] = 0.0;
    }
    tcount = 0;
    
    // q list
    for(i=0; i<num_bins; i++) 
    {
        q_list[i] = qmin + del*i;
        // printf("q = %lf\n", q_list[i]);
    }

    // Grid space 
    idx = 0;
    for (i = 0; i < lx; i++) 
    {
        for (j = 0; j < ly; j++) 
        {
            for (k = 0; k < lz; k++) 
            {
                xyz_grid_center[idx * 3] = i + 0.5;       // Center x
                xyz_grid_center[idx * 3 + 1] = j + 0.5;   // Center y
                xyz_grid_center[idx * 3 + 2] = k + 0.5;   // Center z
                idx++;
            }
        }
    }

    // Each timestep, calculate s(q)
    for(t=tmax-3; t<tmax; t++)
    {
        q_index = 0;

        // Initialize grid sum
        for(i=0; i<lx*ly*lz; i++)
        {
            count_grid[i] = 0;
            count_bin[i] = 0;
        }
        // Save coordinates to grid
        for (i = 0; i < natom; i++) 
        {
            if(type[i] == 1) 
            {
                x_idx = (int)(x[t][i]/1);
                if (x_idx >= 100) x_idx = 0;
                if (x_idx < 0) x_idx = lx - 1;

                y_idx = (int)(y[t][i]/1);
                if (y_idx >= 100) y_idx = 0;
                if (y_idx < 0) y_idx = ly - 1;

                z_idx = (int)(z[t][i]/1);
                if (z_idx >= 20) z_idx = 0;
                if (z_idx < 0) z_idx = lz - 1;

                tot_idx = x_idx * lx*lx + y_idx * ly + z_idx;
                if (tot_idx >= lx*ly*lz) 
                {
                    printf("Error: Out of bounds index %d (x_idx=%d, y_idx=%d, z_idx=%d)\n", 
                        tot_idx, x_idx, y_idx, z_idx);
                } 
                else 
                {
                    count_grid[tot_idx]++;
                    
                }
            }
        }
        // Compute density grid
        for (i = 0; i < lx*ly*lz; i++) 
        {
            if(count_grid[i] > 2) printf("count_grid[%d] = %d\n", i, count_grid[i]);
            dens_grid[i] = count_grid[i] * lx*ly*lz / natom;
            
        }
        
        // Calculate s(q)
        for(qx=0; qx<qbin*10; qx++)
        {
            for(qy=0; qy<qbin*10; qy++)
            {
                for(qz=0; qz<qbin*2; qz++)
                {
                    qxi = exp(qmin+qx*oned_del); qyi = exp(qmin+qy*oned_del); qzi = exp(qmin+qz*oned_del);
                    qmag = sqrt(qxi*qxi + qyi*qyi + qzi*qzi);
                    if(qmag > qmax) continue;
                    else
                    {
                        for(i=0; i<num_bins; i++)
                        {
                            if(qmag > q_list[i] && qmag < q_list[i+1])
                            {
                                count_bin[i]++;
                                q_index = i;
                                break;
                            }
                        }
                    }
                    real_sum = 0.0;
                    imag_sum = 0.0;
                    for(i=0; i<lx*ly*lz; i++)
                    {
                        qdotr[i] = qxi*xyz_grid_center[i*3] + qyi*xyz_grid_center[i*3+1] + qzi*xyz_grid_center[i*3+2];
                        real_sum += cos(qdotr[i]);
                        imag_sum += sin(qdotr[i]);
                    }
                    s_q_value = (real_sum*real_sum + imag_sum*imag_sum);
                    s_q[q_index] += s_q_value;
                    q_real[q_index] += qmag;
                }
            }
            q_index++;
        }
        for(i=0; i<q_index; i++)
        {
            s_q[i] /= count_bin[i];
            q_real[i] /= count_bin[i];
        }
        tcount++;
    }
    // printf("done_sq calc\n");
    for(i=0; i<q_index; i++)
    {
        s_q[i] /= tcount;
        q_real[i] /= tcount;
    }

    // Write average s(q) to file
    FILE *s_q_file = fopen("s_q.txt", "w");
    for(q=0; q<q_index; q++)
    {
        fprintf(s_q_file, "%lf %lf\n", q_real[q], s_q[q]);
        // printf("%lf %lf\n", q_real[q], s_q[q]);
    }
    fclose(s_q_file);
    
    // Free memory
    free(timestep); free(mol); free(type); free(s_q); free(q_list); free(q_real); free(qdotr); free(xyz_grid_center); free(count_grid); free(count_bin); free(dens_grid);
    for(i = 0; i < tmax; i++) free(x[i]);
    free(x);
    for(i = 0; i < tmax; i++) free(y[i]);
    free(y);
    for(i = 0; i < tmax; i++) free(z[i]);
    free(z);
}