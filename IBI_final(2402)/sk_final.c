#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

int main (int argc, const char * argv[])
{
    int i, j, t, q, natom, tmax, qmax, tcount, q_index, dum, id;
    double lx, ly, lz, x1, x2, y1, y2, z1, z2, rx, ry, rz, qmin, qbin, qxi, qyi, qzi, qx, qy, qz;
    double real_sum, imag_sum, r, s_q_value;
    
    // Parameters
    tmax = 22;
    qbin = 2;
    qmin = 0.2;
    qmax = qbin*10*qbin*10*qbin*2;

    FILE *traj = fopen(argv[1], "r");
    fscanf(traj,"ITEM: TIMESTEP\n%d\n", &dum);
    fscanf(traj,"ITEM: NUMBER OF ATOMS\n%d\n", &natom);
    fscanf(traj,"ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
    fscanf(traj,"%lf %lf %lf\n", &x1, &lx, &x2);
    fscanf(traj,"%lf %lf %lf\n", &y1, &ly, &y2);
    fscanf(traj,"%lf %lf %lf\n", &z1, &lz, &z2);
    fscanf(traj,"ITEM: ATOMS id mol type x y z\n");
    fscanf(traj,"ITEM: ATOMS id mol type x y z\n");

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
    double *s_q = calloc(qmax, sizeof(double));
    double *q_list = calloc(qmax, sizeof(double));
    for (i = 0; i < natom; ++i) {
        fscanf(traj, "%d %d %d %lf %lf %lf\n", &id, &mol[i], &type[i], &x[0][i], &y[0][i], &z[0][i]);
    }    
    fclose(traj);
    
    // Read dump file
    traj = fopen(argv[1], "r");
    for(t=0; t<tmax; t++)
    {
        fscanf(traj,"ITEM: TIMESTEP\n%d\n", &dum);
        fscanf(traj,"ITEM: NUMBER OF ATOMS\n%d\n", &natom);
        fscanf(traj,"ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
        fscanf(traj,"%lf %lf %lf\n", &x1, &lx, &x2);
        fscanf(traj,"%lf %lf %lf\n", &y1, &ly, &y2);
        fscanf(traj,"%lf %lf %lf\n", &z1, &lz, &z2);
        fscanf(traj,"ITEM: ATOMS id mol type x y z\n");
        for (i = 0; i < natom; ++i) {
            fscanf(traj, "%d %d %d %lf %lf %lf\n", &id, &mol[i], &type[i], &x[t][i], &y[t][i], &z[t][i]);
        }
        if(x[t][0] == 0) printf("t=%d\n", t);
    }
    fclose(traj);
   // Initialize s(q)
    for(q=0; q<qmax; q++)
    {
        s_q[q] = 0.0;
    }
    tcount = 0;
    //calculate s(q)
    // q list
    q_index = 0;
    for(qx=-qbin*5; qx<qbin*5; qx++)
    {
        for(qy=-qbin*5; qy<qbin*5; qy++)
        {
            for(qz=0; qz<qbin*2; qz++)
            // for(qz = -qbin*5; qz < qbin*5; qz++)
            {
                qxi = qmin*qx; qyi = qmin*qy; qzi =  qmin*qz;
                q_list[q_index] = sqrt(qxi*qxi + qyi*qyi + qzi*qzi);
                q_index++;
            }
        }
    }
    for(t=(int)(tmax*0); t<tmax; t++)
    {
        q_index = 0;
        for(qx=-qbin*5; qx<qbin*5; qx++)
        {
            for(qy=-qbin*5; qy<qbin*5; qy++)
            {
                for(qz=0; qz<qbin*2; qz++)
                // for(qz = -qbin*5; qz < qbin*5; qz++)
                {
                    qxi = qmin*qx; qyi = qmin*qy; qzi = qmin*qz;
                    real_sum = 0.0;
                    imag_sum = 0.0;
                    for(i=0; i<natom; i++)
                    {
                        if(type[i] == 2)
                        {
                            rx = (x[t][i]) * qxi;
                            ry = (y[t][i]) * qyi;
                            rz = (z[t][i]) * qzi;
                            r = (rx*rx + ry*ry + rz*rz);
                            real_sum += cos(r);
                            imag_sum += sin(r);
                        }
                    }
                    s_q_value = ((real_sum*real_sum + imag_sum*imag_sum));
                    s_q[q_index] += s_q_value/natom;
                    q_index++;
                }
            }
        }
        tcount++;

        FILE *s_q_file = fopen("s_q.txt", "a");
        fprintf(s_q_file, "time = %d\n", tcount);
        for(q=0; q<q_index; q++)
        {
            fprintf(s_q_file, "%lf %lf\n", q_list[q], s_q[q]);
            // printf("%lf %lf\n", q_list[q], s_q[q]);
        }
        fclose(s_q_file);
    }
    // printf("done_sq calc\n");
    for(q=0; q<q_index; q++)
    {
        s_q[q] /= tcount;
    }

    // Write average s(q) to file
    FILE *s_q_file = fopen("s_q.txt", "w");
    for(q=0; q<q_index; q++)
    {
        fprintf(s_q_file, "%lf %lf\n", q_list[q], s_q[q]);
        printf("%lf %lf\n", q_list[q], s_q[q]);
    }
    fclose(s_q_file);
    
    
    // Free memory
    free(timestep); free(mol); free(type); free(s_q); free(q_list);
    for(i = 0; i < tmax; i++) free(x[i]);
    free(x);
    for(i = 0; i < tmax; i++) free(y[i]);
    free(y);
    for(i = 0; i < tmax; i++) free(z[i]);
    free(z);
}