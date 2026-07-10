#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

int main (int argc, const char * argv[])
{
    
    int natom_total, dint, tmax;
    double x, y, z;
    FILE *traj = fopen("traj_blend.xyz", "r");

    FILE *coord = fopen("coord.xyz", "w");
    natom_total = atoi(argv[1]);
    tmax = atoi(argv[2]);

    for(int t=0; t<tmax; t++)
    {
        fscanf(traj, "%d\n", &dint); // Read the number of atoms in the last configuration
        fscanf(traj, "%d 1.500000\n", &dint); // Skip the line containing cell information
        for (int i = 0; i < natom_total; ++i) {
            fscanf(traj, "%d %lf %lf %lf\n", &dint, &x, &y, &z);

            if(t==(tmax-1)) fprintf(coord, "%lf %lf %lf\n", x, y, z); // Write atom information to coord.xyz
        }
    }
    fclose(coord);
    fclose(traj);

}
