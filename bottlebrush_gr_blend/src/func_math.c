#include "header.h"

// get the mod between two integers (guarantee only positive number result)
int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}

// calculate the position (scaler) due to periodic boundary condition
double periodic(double ri, double bl)
{
    if (ri < 0)
        ri += bl;
    else if (ri > bl)
        ri -= bl;

    return ri;
}

// distance calculation that accounts for periodic boundary condition
double dist_calc(double ra_x, double rb_x, double ra_y, double rb_y, double ra_z, double rb_z, double bl)
{
    double dpos[3], bl_half;

    bl_half = bl / 2.0;

    dpos[0] = ra_x - rb_x;
    dpos[1] = ra_y - rb_y;
    dpos[2] = ra_z - rb_z;

    if (fabs(dpos[0]) > bl_half)
        dpos[0] = fabs(dpos[0]) - bl;
    if (fabs(dpos[1]) > bl_half)
        dpos[1] = fabs(dpos[1]) - bl;
    if (fabs(dpos[2]) > bl_half)
        dpos[2] = fabs(dpos[2]) - bl;

    return sqrt(dpos[0] * dpos[0] + dpos[1] * dpos[1] + dpos[2] * dpos[2]);
}

// distance (squared) that accounts for periodic boundary condition
double dist_sq_calc(double ra_x, double rb_x, double ra_y, double rb_y, double ra_z, double rb_z, double bl)
{
    double dpos[3], bl_half;

    bl_half = bl / 2.0;

    dpos[0] = ra_x - rb_x;
    dpos[1] = ra_y - rb_y;
    dpos[2] = ra_z - rb_z;

    if (fabs(dpos[0]) > bl_half)
        dpos[0] = fabs(dpos[0]) - bl;
    if (fabs(dpos[1]) > bl_half)
        dpos[1] = fabs(dpos[1]) - bl;
    if (fabs(dpos[2]) > bl_half)
        dpos[2] = fabs(dpos[2]) - bl;

    return dpos[0] * dpos[0] + dpos[1] * dpos[1] + dpos[2] * dpos[2];
}

// get the periodic image of a position (1d)
double simbox_1dgetimage(double rx, double bl)
{
    double res = 0;
    res = rx - bl * round(rx / bl);
    return res;
}

// metropolis criteria determination of monte carlo moves
int metro_crit(double enrg_diff, long *idum)
{
    int res = 0;
    double ran, prob;
    if (enrg_diff <= 0)
    {
        res = 1;
    }
    else if ((enrg_diff > 0.0) && (enrg_diff <= 40.0))
    {
        ran = ran1(idum);
        prob = exp(-enrg_diff);
        if (ran <= prob)
        {
            res = 1;
        }
    }

    return res;
}

// gr_prep
void gr_prep(int type, int NISC1, int NISC2, int tot_isc1, int tot_isc2, double **com1, double **com2, double box_size, int nbin, int *histogram)
{
    int i, j, bin;
    double dx, dy, dz, dr, r, rx1, ry1, rz1, rx2, ry2, rz2, res;
    res = 0;
    if(type == 0) // Same type, com-com, com_mid-com_mid pairs
    {
        for (i = 0; i < tot_isc1 - 1; i++) 
        {
            for (j = i + 1; j < tot_isc2; j++) 
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
                if (bin < nbin) histogram[bin]++;
            }
        }
    }
    else if (type ==1) //different type, com-com, com_mid-com_mid, com-com_mid pairs 
    {
        for (i = 0; i < tot_isc1; i++) 
        {
            for (j = 0; j < tot_isc2; j++) 
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
                if (bin < nbin) histogram[bin]++;
            }
        }
    }
    else if (type == 2) // same type, com-com_mid pairs
    {
        for (i = 0; i < tot_isc1; i++) 
        {
            for (j = 0; j < tot_isc2; j++) 
            {
                if((int)i/(NISC1) == (int)j/(NISC2) && (abs((int)i%(NISC1)-(int)j%(NISC2)) < 2)) continue;
                else
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
                    if (bin < nbin) histogram[bin]++;
                }
            }
        }
    }
}