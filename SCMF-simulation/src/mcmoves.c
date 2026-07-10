#include "header.h"
// monte carlo move for a single atom displacement
// UPDATE: bond energy should be divided by 2 species
void mc_atom_displ(int iatom, double displ, int pm_type, unsigned long long *naccept, double ****density_grids, double density,
                   int *particle_type, double **coords, int **atom_bond_list, double box_size, double grid_size, double rho_norm,
                   int maxsite_1d, double chi, double kappa, double kspring, long long int natom_total,
                   double grid_shift[3], long *idum)
{
    if (pm_type == 0)
    {
        double rx, ry, rz;
        int xgrid, ygrid, zgrid;
        int atom_type;
        double kappaN, pre_factor;

        atom_type = particle_type[iatom];
        rx = coords[iatom][0];
        ry = coords[iatom][1];
        rz = coords[iatom][2];
        // identify the grid that the particle is originally in, taking a grid shift into account
        xgrid = rx / grid_size + grid_shift[0];
        xgrid = xgrid % maxsite_1d;
        ygrid = ry / grid_size + grid_shift[1];
        ygrid = ygrid % maxsite_1d;
        zgrid = rz / grid_size + grid_shift[2];
        zgrid = zgrid % maxsite_1d;

        double theta, phi;
        double rx_new, ry_new, rz_new;
        int xgrid_new, ygrid_new, zgrid_new;
        // Randomly get theta angle(xy plane)
        theta = ran1(idum) * 2.0 * M_PI;
        // randomly get cos() value ranging from -1 to 1 to set 2Pi rotation angle from z axis
        phi = acos(2.0 * ran1(idum) - 1.0);
        // move the particle by a vector with length displ
        rx_new = rx + displ * cos(theta) * sin(phi);
        ry_new = ry + displ * sin(theta) * sin(phi);
        rz_new = rz + displ * cos(phi);
        // wrap around the periodic boundary condition
        rx_new = periodic(rx_new, box_size);
        ry_new = periodic(ry_new, box_size);
        rz_new = periodic(rz_new, box_size);
        // identify the grid that the particle is moved to, taking a grid shift into account
        xgrid_new = rx_new / grid_size + grid_shift[0];
        xgrid_new = xgrid_new % maxsite_1d;
        ygrid_new = ry_new / grid_size + grid_shift[1];
        ygrid_new = ygrid_new % maxsite_1d;
        zgrid_new = rz_new / grid_size + grid_shift[2];
        zgrid_new = zgrid_new % maxsite_1d;
        // printf("x y z: %d %d %d\n", xgrid, ygrid, zgrid);
        // printf("x_new y_new z_new: %d %d %d\n", xgrid_new, ygrid_new, zgrid_new);

        double phiA_new, phiB_new;
        double phiA_old_from, phiB_old_from, phiA_old_to, phiB_old_to;
        double phiA_new_from, phiB_new_from, phiA_new_to, phiB_new_to;
        double enrg_local_old, enrg_local_new;
        double enrg_bond_old, enrg_bond_new;
        double enrg_diff;
        int iaccept;

        // Update the new density, ONLY WHEN THE GRID CHANGES!
        if (xgrid_new != xgrid || ygrid_new != xgrid || zgrid_new != zgrid)
        {
            // record the old density of the grid that this particle is:
            // moving from (_from) and moving to (_to)
            phiA_old_from = density_grids[0][xgrid][ygrid][zgrid];
            phiB_old_from = density_grids[1][xgrid][ygrid][zgrid];
            phiA_old_to = density_grids[0][xgrid_new][ygrid_new][zgrid_new];
            phiB_old_to = density_grids[1][xgrid_new][ygrid_new][zgrid_new];

            enrg_local_old = calc_singlemesh_energy(phiA_old_from, phiB_old_from, chi, kappa, density) +
                             calc_singlemesh_energy(phiA_old_to, phiB_old_to, chi, kappa, density);
            if (atom_type == 0)
            {
                phiA_new_from = phiA_old_from - 1 / rho_norm;
                phiA_new_to = phiA_old_to + 1 / rho_norm;
                phiB_new_from = phiB_old_from;
                phiB_new_to = phiB_old_to;
            }
            else
            {
                phiB_new_from = phiB_old_from - 1 / rho_norm;
                phiB_new_to = phiB_old_to + 1 / rho_norm;
                phiA_new_from = phiA_old_from;
                phiA_new_to = phiA_old_to;
            }
            // calculate the single mesh energy
            enrg_local_new = calc_singlemesh_energy(phiA_new_from, phiB_new_from, chi, kappa, density) +
                             calc_singlemesh_energy(phiA_new_to, phiB_new_to, chi, kappa, density);
            // Calculate bond energy
            enrg_bond_old = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);

            // move the particle to the new position and calculate the new bond energy
            coords[iatom][0] = rx_new;
            coords[iatom][1] = ry_new;
            coords[iatom][2] = rz_new;
            enrg_bond_new = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);

            enrg_diff = (enrg_local_new - enrg_local_old) +
                        (enrg_bond_new - enrg_bond_old);

            // Use metropolis criteria to determine whether to accept the move
            iaccept = metro_crit(enrg_diff, idum);
            // For debug:
            // printf("phiA_old_from = %lf, phiA_old_to = %lf\n", phiA_old_from, phiA_old_to);
            // printf("phiB_old_from = %lf, phiB_old_to = %lf\n", phiB_old_from, phiB_old_to);
            // printf("enrg_local_old = %lf, enrg_local_new = %lf\n, \t enrg_bond_old = %lf, enrg_bond_new = %lf\n, \t enrg_diff = %lf, accept = %d\n", enrg_local_old, enrg_local_new, enrg_bond_old, enrg_bond_new, enrg_diff, iaccept);
            *naccept += iaccept;
            if (iaccept == 0)
            {
                // move back the particle if not accepted
                coords[iatom][0] = rx;
                coords[iatom][1] = ry;
                coords[iatom][2] = rz;
            }
            else
            {
                // update the density profile if accepted
                if (atom_type == 0)
                {
                    density_grids[0][xgrid][ygrid][zgrid] = phiA_new_from;
                    density_grids[0][xgrid_new][ygrid_new][zgrid_new] = phiA_new_to;
                }
                else
                {
                    density_grids[1][xgrid][ygrid][zgrid] = phiB_new_from;
                    density_grids[1][xgrid_new][ygrid_new][zgrid_new] = phiB_new_to;
                }
            }
        }
        // only calculate bond energy if the grid is not changed
        else
        {
            // Calculate bond energy
            enrg_bond_old = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);

            // move the particle to the new position and calculate the new bond energy
            coords[iatom][0] = rx_new;
            coords[iatom][1] = ry_new;
            coords[iatom][2] = rz_new;
            enrg_bond_new = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);
            enrg_diff = enrg_bond_new - enrg_bond_old;

            // Use metropolis criteria to determine whether to accept the move
            iaccept = metro_crit(enrg_diff, idum);
            *naccept += iaccept;
            if (iaccept == 0)
            {
                // move back the particle if not accepted
                coords[iatom][0] = rx;
                coords[iatom][1] = ry;
                coords[iatom][2] = rz;
            }
        }
    }

    if (pm_type == 1)
    {
        double rx, ry, rz;
        int xgrid, ygrid, zgrid;
        int atom_type;
        double kappaN, pre_factor;        
        double xc, yc, zc;
        double dx, dy, dz;
        double tx, ty, tz;
        double dx2, dy2, dz2;
        double tx2, ty2, tz2;
        int xneigh, yneigh, zneigh;

        atom_type = particle_type[iatom];
        rx = coords[iatom][0];
        ry = coords[iatom][1];
        rz = coords[iatom][2];

        // identify the grid that the particle is in and the corresponding grid center;
        // taking a grid shift into account
        xgrid = rx / grid_size + grid_shift[0];
        xc = ((double)xgrid + 0.5 - grid_shift[0]) * grid_size;
        xgrid = mod(xgrid, maxsite_1d);
        ygrid = ry / grid_size + grid_shift[1];
        yc = ((double)ygrid + 0.5 - grid_shift[1]) * grid_size;
        ygrid = mod(ygrid, maxsite_1d);
        zgrid = rz / grid_size + grid_shift[2];
        zc = ((double)zgrid + 0.5 - grid_shift[2]) * grid_size;
        zgrid = mod(zgrid, maxsite_1d);

        dx = (rx - xc) / grid_size;
        dy = (ry - yc) / grid_size;
        dz = (rz - zc) / grid_size;

        tx = 1 - fabs(dx);
        ty = 1 - fabs(dy);
        tz = 1 - fabs(dz);

        xneigh = xgrid + ((dx > 0) ? 1 : -1);
        xneigh = mod(xneigh, maxsite_1d);
        yneigh = ygrid + ((dy > 0) ? 1 : -1);
        yneigh = mod(yneigh, maxsite_1d);
        zneigh = zgrid + ((dz > 0) ? 1 : -1);
        zneigh = mod(zneigh, maxsite_1d);
        //if(xgrid < 0 || ygrid < 0 || zgrid < 0 || xneigh < 0 || yneigh < 0 || zneigh < 0)
        //{
        //    printf("xgrid = %d, ygrid = %d, zgrid = %d\n", xgrid, ygrid, zgrid);
        //    printf("xneigh = %d, yneigh = %d, zneigh = %d\n", xneigh, yneigh, zneigh);
        //}
        //if(dx <0 || dy < 0 || dz < 0 || tx < 0 || ty < 0 || tz < 0)
        //{
        //    printf("dx = %f, dy = %f, dz = %f\n", dx, dy, dz);
        //    printf("tx = %f, ty = %f, tz = %f\n", tx, ty, tz);
        //}
        // Randomly displace the particle
        double theta, phi;
        double rx_new, ry_new, rz_new;
        int xgrid_new, ygrid_new, zgrid_new, xneigh_new, yneigh_new, zneigh_new;

        theta = ran1(idum) * 2.0 * M_PI;
        // randomly get cos() value ranging from -1 to 1 to set 2Pi rotation angle from z axis
        phi = acos(2.0 * ran1(idum) - 1.0);
        // move the particle by a vector with length displ
        rx_new = rx + displ * cos(theta) * sin(phi);
        ry_new = ry + displ * sin(theta) * sin(phi);
        rz_new = rz + displ * cos(phi);
        // wrap around the periodic boundary condition
        rx_new = periodic(rx_new, box_size);
        ry_new = periodic(ry_new, box_size);
        rz_new = periodic(rz_new, box_size);

        // identify the grid that the particle is in and the corresponding grid center;
        xgrid_new = rx_new / grid_size + grid_shift[0];
        xc = ((double)xgrid_new + 0.5 - grid_shift[0]) * grid_size;
        xgrid_new = mod(xgrid_new, maxsite_1d);
        ygrid_new = ry_new / grid_size + grid_shift[1];
        yc = ((double)ygrid_new + 0.5 - grid_shift[1]) * grid_size;
        ygrid_new = mod(ygrid_new, maxsite_1d);
        zgrid_new = rz_new / grid_size + grid_shift[2];
        zc = ((double)zgrid_new + 0.5 - grid_shift[2]) * grid_size;
        zgrid_new = mod(zgrid_new, maxsite_1d);

        dx2 = (rx_new - xc) / grid_size;
        dy2 = (ry_new - yc) / grid_size;
        dz2 = (rz_new - zc) / grid_size;

        tx2 = 1 - fabs(dx2);
        ty2 = 1 - fabs(dy2);
        tz2 = 1 - fabs(dz2);

        xneigh_new = xgrid_new + ((dx2 > 0) ? 1 : -1);
        xneigh_new = mod(xneigh_new, maxsite_1d);
        yneigh_new = ygrid_new + ((dy2 > 0) ? 1 : -1);
        yneigh_new = mod(yneigh_new, maxsite_1d);
        zneigh_new = zgrid_new + ((dz2 > 0) ? 1 : -1);
        zneigh_new = mod(zneigh_new, maxsite_1d);

        double phiA_new, phiB_new;
        double phiA_old_from, phiB_old_from, phiA_old_to, phiB_old_to;
        double phiA_new_from, phiB_new_from, phiA_new_to, phiB_new_to;
        double enrg_local_old, enrg_local_new;
        double enrg_bond_old, enrg_bond_new;
        double enrg_diff;
        int iaccept;
        double phiA_old_from1, phiA_old_from2, phiA_old_from3, phiA_old_from4;
        double phiA_old_from5, phiA_old_from6, phiA_old_from7, phiA_old_from8;
        double phiA_old_to1, phiA_old_to2, phiA_old_to3, phiA_old_to4;
        double phiA_old_to5, phiA_old_to6, phiA_old_to7, phiA_old_to8;
        double phiB_old_from1, phiB_old_from2, phiB_old_from3, phiB_old_from4;
        double phiB_old_from5, phiB_old_from6, phiB_old_from7, phiB_old_from8;
        double phiB_old_to1, phiB_old_to2, phiB_old_to3, phiB_old_to4;
        double phiB_old_to5, phiB_old_to6, phiB_old_to7, phiB_old_to8;
        double phiA_new_from1, phiA_new_from2, phiA_new_from3, phiA_new_from4;
        double phiA_new_from5, phiA_new_from6, phiA_new_from7, phiA_new_from8;
        double phiA_new_to1, phiA_new_to2, phiA_new_to3, phiA_new_to4;
        double phiA_new_to5, phiA_new_to6, phiA_new_to7, phiA_new_to8;
        double phiB_new_from1, phiB_new_from2, phiB_new_from3, phiB_new_from4;
        double phiB_new_from5, phiB_new_from6, phiB_new_from7, phiB_new_from8;
        double phiB_new_to1, phiB_new_to2, phiB_new_to3, phiB_new_to4;
        double phiB_new_to5, phiB_new_to6, phiB_new_to7, phiB_new_to8;

        //Update new density when the grid changes
        if(xgrid_new == xgrid && ygrid_new == ygrid && zgrid_new == zgrid && xneigh == xneigh_new && yneigh == yneigh_new && zneigh == zneigh_new)
        {
            //calculate bond energy
            enrg_bond_old = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);
            coords[iatom][0] = rx_new;
            coords[iatom][1] = ry_new;
            coords[iatom][2] = rz_new;
            enrg_bond_new = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);
            //calculate the energy difference
            enrg_diff = enrg_bond_new - enrg_bond_old;        
            // Use metroplis algorithm to decide whether to accept the move
            iaccept = metro_crit(enrg_diff, idum);
            *naccept += iaccept;
            if(iaccept == 0)
            {
                coords[iatom][0] = rx;
                coords[iatom][1] = ry;
                coords[iatom][2] = rz;
            }
        }
        else
        {
            phiA_old_from1 = density_grids[0][xgrid][ygrid][zgrid];
            phiA_old_from2 = density_grids[0][xneigh][ygrid][zgrid];
            phiA_old_from3 = density_grids[0][xgrid][yneigh][zgrid];
            phiA_old_from4 = density_grids[0][xgrid][ygrid][zneigh];
            phiA_old_from5 = density_grids[0][xneigh][yneigh][zgrid];
            phiA_old_from6 = density_grids[0][xneigh][ygrid][zneigh];
            phiA_old_from7 = density_grids[0][xgrid][yneigh][zneigh];
            phiA_old_from8 = density_grids[0][xneigh][yneigh][zneigh];

            phiA_old_to1 = density_grids[0][xgrid_new][ygrid_new][zgrid_new];
            phiA_old_to2 = density_grids[0][xneigh_new][ygrid_new][zgrid_new];
            phiA_old_to3 = density_grids[0][xgrid_new][yneigh_new][zgrid_new];
            phiA_old_to4 = density_grids[0][xgrid_new][ygrid_new][zneigh_new];
            phiA_old_to5 = density_grids[0][xneigh_new][yneigh_new][zgrid_new];
            phiA_old_to6 = density_grids[0][xneigh_new][ygrid_new][zneigh_new];
            phiA_old_to7 = density_grids[0][xgrid_new][yneigh_new][zneigh_new];
            phiA_old_to8 = density_grids[0][xneigh_new][yneigh_new][zneigh_new];

            phiB_old_from1 = density_grids[1][xgrid][ygrid][zgrid];
            phiB_old_from2 = density_grids[1][xneigh][ygrid][zgrid];
            phiB_old_from3 = density_grids[1][xgrid][yneigh][zgrid];
            phiB_old_from4 = density_grids[1][xgrid][ygrid][zneigh];
            phiB_old_from5 = density_grids[1][xneigh][yneigh][zgrid];
            phiB_old_from6 = density_grids[1][xneigh][ygrid][zneigh];
            phiB_old_from7 = density_grids[1][xgrid][yneigh][zneigh];
            phiB_old_from8 = density_grids[1][xneigh][yneigh][zneigh];

            phiB_old_to1 = density_grids[1][xgrid_new][ygrid_new][zgrid_new];
            phiB_old_to2 = density_grids[1][xneigh_new][ygrid_new][zgrid_new];
            phiB_old_to3 = density_grids[1][xgrid_new][yneigh_new][zgrid_new];
            phiB_old_to4 = density_grids[1][xgrid_new][ygrid_new][zneigh_new];
            phiB_old_to5 = density_grids[1][xneigh_new][yneigh_new][zgrid_new];
            phiB_old_to6 = density_grids[1][xneigh_new][ygrid_new][zneigh_new];
            phiB_old_to7 = density_grids[1][xgrid_new][yneigh_new][zneigh_new];
            phiB_old_to8 = density_grids[1][xneigh_new][yneigh_new][zneigh_new];

            enrg_local_old = calc_singlemesh_energy(phiA_old_from1, phiB_old_from1, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from2, phiB_old_from2, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from3, phiB_old_from3, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from4, phiB_old_from4, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from5, phiB_old_from5, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from6, phiB_old_from6, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from7, phiB_old_from7, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from8, phiB_old_from8, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to1, phiB_old_to1, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to2, phiB_old_to2, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to3, phiB_old_to3, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to4, phiB_old_to4, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to5, phiB_old_to5, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to6, phiB_old_to6, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to7, phiB_old_to7, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to8, phiB_old_to8, chi, kappa, density);

            //Here must have bug...
            if (atom_type == 0)
            {
                phiA_new_from1 = phiA_old_from1 - tx*ty*tz/rho_norm;
                phiA_new_from2 = phiA_old_from2 - (1-tx)*ty*tz/rho_norm;
                phiA_new_from3 = phiA_old_from3 - tx*(1-ty)*tz/rho_norm;
                phiA_new_from4 = phiA_old_from4 - tx*ty*(1-tz)/rho_norm;
                phiA_new_from5 = phiA_old_from5 - (1-tx)*(1-ty)*tz/rho_norm;
                phiA_new_from6 = phiA_old_from6 - (1-tx)*ty*(1-tz)/rho_norm;
                phiA_new_from7 = phiA_old_from7 - tx*(1-ty)*(1-tz)/rho_norm;
                phiA_new_from8 = phiA_old_from8 - (1-tx)*(1-ty)*(1-tz)/rho_norm;
                
                phiA_new_to1 = phiA_old_to1 + tx*ty*tz/rho_norm;
                phiA_new_to2 = phiA_old_to2 + (1-tx)*ty*tz/rho_norm;
                phiA_new_to3 = phiA_old_to3 + tx*(1-ty)*tz/rho_norm;
                phiA_new_to4 = phiA_old_to4 + tx*ty*(1-tz)/rho_norm;
                phiA_new_to5 = phiA_old_to5 + (1-tx)*(1-ty)*tz/rho_norm;
                phiA_new_to6 = phiA_old_to6 + (1-tx)*ty*(1-tz)/rho_norm;
                phiA_new_to7 = phiA_old_to7 + tx*(1-ty)*(1-tz)/rho_norm;
                phiA_new_to8 = phiA_old_to8 + (1-tx)*(1-ty)*(1-tz)/rho_norm;

                phiB_new_from1 = phiB_old_from1;
                phiB_new_from2 = phiB_old_from2;
                phiB_new_from3 = phiB_old_from3;
                phiB_new_from4 = phiB_old_from4;
                phiB_new_from5 = phiB_old_from5;
                phiB_new_from6 = phiB_old_from6;
                phiB_new_from7 = phiB_old_from7;
                phiB_new_from8 = phiB_old_from8;

                phiB_new_to1 = phiB_old_to1;
                phiB_new_to2 = phiB_old_to2;
                phiB_new_to3 = phiB_old_to3;
                phiB_new_to4 = phiB_old_to4;
                phiB_new_to5 = phiB_old_to5;
                phiB_new_to6 = phiB_old_to6;
                phiB_new_to7 = phiB_old_to7;
                phiB_new_to8 = phiB_old_to8;
            }
            else
            {
                phiB_new_from1 = phiB_old_from1 - tx*ty*tz/rho_norm;
                phiB_new_from2 = phiB_old_from2 - (1-tx)*ty*tz/rho_norm;
                phiB_new_from3 = phiB_old_from3 - tx*(1-ty)*tz/rho_norm;
                phiB_new_from4 = phiB_old_from4 - tx*ty*(1-tz)/rho_norm;
                phiB_new_from5 = phiB_old_from5 - (1-tx)*(1-ty)*tz/rho_norm;
                phiB_new_from6 = phiB_old_from6 - (1-tx)*ty*(1-tz)/rho_norm;
                phiB_new_from7 = phiB_old_from7 - tx*(1-ty)*(1-tz)/rho_norm;
                phiB_new_from8 = phiB_old_from8 - (1-tx)*(1-ty)*(1-tz)/rho_norm;

                phiB_new_to1 = phiB_old_to1 + tx*ty*tz/rho_norm;
                phiB_new_to2 = phiB_old_to2 + (1-tx)*ty*tz/rho_norm;
                phiB_new_to3 = phiB_old_to3 + tx*(1-ty)*tz/rho_norm;
                phiB_new_to4 = phiB_old_to4 + tx*ty*(1-tz)/rho_norm;
                phiB_new_to5 = phiB_old_to5 + (1-tx)*(1-ty)*tz/rho_norm;
                phiB_new_to6 = phiB_old_to6 + (1-tx)*ty*(1-tz)/rho_norm;
                phiB_new_to7 = phiB_old_to7 + tx*(1-ty)*(1-tz)/rho_norm;
                phiB_new_to8 = phiB_old_to8 + (1-tx)*(1-ty)*(1-tz)/rho_norm;

                phiA_new_from1 = phiA_old_from1;
                phiA_new_from2 = phiA_old_from2;
                phiA_new_from3 = phiA_old_from3;
                phiA_new_from4 = phiA_old_from4;
                phiA_new_from5 = phiA_old_from5;
                phiA_new_from6 = phiA_old_from6;
                phiA_new_from7 = phiA_old_from7;
                phiA_new_from8 = phiA_old_from8;

                phiA_new_to1 = phiA_old_to1;
                phiA_new_to2 = phiA_old_to2;
                phiA_new_to3 = phiA_old_to3;
                phiA_new_to4 = phiA_old_to4;
                phiA_new_to5 = phiA_old_to5;
                phiA_new_to6 = phiA_old_to6;
                phiA_new_to7 = phiA_old_to7;
                phiA_new_to8 = phiA_old_to8;
            }
            // calculate the single mesh energy
            enrg_local_new = calc_singlemesh_energy(phiA_new_from1, phiB_new_from1, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from2, phiB_new_from2, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from3, phiB_new_from3, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from4, phiB_new_from4, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from5, phiB_new_from5, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from6, phiB_new_from6, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from7, phiB_new_from7, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from8, phiB_new_from8, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to1, phiB_new_to1, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to2, phiB_new_to2, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to3, phiB_new_to3, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to4, phiB_new_to4, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to5, phiB_new_to5, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to6, phiB_new_to6, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to7, phiB_new_to7, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to8, phiB_new_to8, chi, kappa, density);
            //printf("%f %f %f %f %f %f %f %f\n",calc_singlemesh_energy(phiA_new_from1, phiB_new_from1, chi, kappa, density), calc_singlemesh_energy(phiA_new_from2, phiB_new_from2, chi, kappa, density), calc_singlemesh_energy(phiA_new_from3, phiB_new_from3, chi, kappa, density), calc_singlemesh_energy(phiA_new_from4, phiB_new_from4, chi, kappa, density), calc_singlemesh_energy(phiA_new_from5, phiB_new_from5, chi, kappa, density), calc_singlemesh_energy(phiA_new_from6, phiB_new_from6, chi, kappa, density), calc_singlemesh_energy(phiA_new_from7, phiB_new_from7, chi, kappa, density), calc_singlemesh_energy(phiA_new_from8, phiB_new_from8, chi, kappa, density));
            //calculate bond energy
            enrg_bond_old = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);
            coords[iatom][0] = rx_new;
            coords[iatom][1] = ry_new;
            coords[iatom][2] = rz_new;
            enrg_bond_new = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);
            //calculate the energy difference
            enrg_diff = (enrg_local_new - enrg_local_old) + (enrg_bond_new - enrg_bond_old);  
            //printf("enrg_diff = %f\n", enrg_diff);      
            // Use metroplis algorithm to decide whether to accept the move
            iaccept = metro_crit(enrg_diff, idum);
            //printf("iaccept = %d\n", iaccept);
            *naccept += iaccept;
            if(iaccept == 0)
            {
                coords[iatom][0] = rx;
                coords[iatom][1] = ry;
                coords[iatom][2] = rz;
            }
            else
            {
                coords[iatom][0] = rx_new;
                coords[iatom][1] = ry_new;
                coords[iatom][2] = rz_new;
                
                if(atom_type == 0)
                {
                    density_grids[0][xgrid][ygrid][zgrid] = phiA_new_from1;
                    density_grids[0][xneigh][ygrid][zgrid] = phiA_new_from2;
                    density_grids[0][xgrid][yneigh][zgrid] = phiA_new_from3;
                    density_grids[0][xgrid][ygrid][zneigh] = phiA_new_from4;
                    density_grids[0][xneigh][yneigh][zgrid] = phiA_new_from5;
                    density_grids[0][xneigh][ygrid][zneigh] = phiA_new_from6;
                    density_grids[0][xgrid][yneigh][zneigh] = phiA_new_from7;
                    density_grids[0][xneigh][yneigh][zneigh] = phiA_new_from8;

                    density_grids[0][xgrid_new][ygrid_new][zgrid_new] = phiA_new_to1;
                    density_grids[0][xneigh_new][ygrid_new][zgrid_new] = phiA_new_to2;
                    density_grids[0][xgrid_new][yneigh_new][zgrid_new] = phiA_new_to3;
                    density_grids[0][xgrid_new][ygrid_new][zneigh_new] = phiA_new_to4;
                    density_grids[0][xneigh_new][yneigh_new][zgrid_new] = phiA_new_to5;
                    density_grids[0][xneigh_new][ygrid_new][zneigh_new] = phiA_new_to6;
                    density_grids[0][xgrid_new][yneigh_new][zneigh_new] = phiA_new_to7;
                    density_grids[0][xneigh_new][yneigh_new][zneigh_new] = phiA_new_to8;
                }
                else
                {
                    density_grids[1][xgrid][ygrid][zgrid] = phiB_new_from1;
                    density_grids[1][xneigh][ygrid][zgrid] = phiB_new_from2;
                    density_grids[1][xgrid][yneigh][zgrid] = phiB_new_from3;
                    density_grids[1][xgrid][ygrid][zneigh] = phiB_new_from4;
                    density_grids[1][xneigh][yneigh][zgrid] = phiB_new_from5;
                    density_grids[1][xneigh][ygrid][zneigh] = phiB_new_from6;
                    density_grids[1][xgrid][yneigh][zneigh] = phiB_new_from7;
                    density_grids[1][xneigh][yneigh][zneigh] = phiB_new_from8;

                    density_grids[1][xgrid_new][ygrid_new][zgrid_new] = phiB_new_to1;
                    density_grids[1][xneigh_new][ygrid_new][zgrid_new] = phiB_new_to2;
                    density_grids[1][xgrid_new][yneigh_new][zgrid_new] = phiB_new_to3;
                    density_grids[1][xgrid_new][ygrid_new][zneigh_new] = phiB_new_to4;
                    density_grids[1][xneigh_new][yneigh_new][zgrid_new] = phiB_new_to5;
                    density_grids[1][xneigh_new][ygrid_new][zneigh_new] = phiB_new_to6;
                    density_grids[1][xgrid_new][yneigh_new][zneigh_new] = phiB_new_to7;
                    density_grids[1][xneigh_new][yneigh_new][zneigh_new] = phiB_new_to8;

                }
            
            }
        }

    }
}

void mc_atom_displ_bond(int iatom, double displ, int pm_type, unsigned long long *naccept, double ****density_grids, double density,
                   int *particle_type, double **coords, int **atom_bond_list, double box_size, double grid_size, double rho_norm,
                   int maxsite_1d, double chi, double kappa, double kspring, long long int natom_total,
                   double grid_shift[3], long *idum)
{
    if (pm_type == 0)
    {
        double rx, ry, rz;
        int xgrid, ygrid, zgrid;
        int atom_type;
        double kappaN, pre_factor;

        atom_type = particle_type[iatom];
        rx = coords[iatom][0];
        ry = coords[iatom][1];
        rz = coords[iatom][2];
        // identify the grid that the particle is originally in, taking a grid shift into account
        xgrid = rx / grid_size + grid_shift[0];
        xgrid = xgrid % maxsite_1d;
        ygrid = ry / grid_size + grid_shift[1];
        ygrid = ygrid % maxsite_1d;
        zgrid = rz / grid_size + grid_shift[2];
        zgrid = zgrid % maxsite_1d;

        double theta, phi;
        double rx_new, ry_new, rz_new;
        int xgrid_new, ygrid_new, zgrid_new;
        // Randomly get theta angle(xy plane)
        theta = ran1(idum) * 2.0 * M_PI;
        // randomly get cos() value ranging from -1 to 1 to set 2Pi rotation angle from z axis
        phi = acos(2.0 * ran1(idum) - 1.0);
        // move the particle by a vector with length displ
        rx_new = rx + displ * cos(theta) * sin(phi);
        ry_new = ry + displ * sin(theta) * sin(phi);
        rz_new = rz + displ * cos(phi);
        // wrap around the periodic boundary condition
        rx_new = periodic(rx_new, box_size);
        ry_new = periodic(ry_new, box_size);
        rz_new = periodic(rz_new, box_size);
        // identify the grid that the particle is moved to, taking a grid shift into account
        xgrid_new = rx_new / grid_size + grid_shift[0];
        xgrid_new = xgrid_new % maxsite_1d;
        ygrid_new = ry_new / grid_size + grid_shift[1];
        ygrid_new = ygrid_new % maxsite_1d;
        zgrid_new = rz_new / grid_size + grid_shift[2];
        zgrid_new = zgrid_new % maxsite_1d;
        // printf("x y z: %d %d %d\n", xgrid, ygrid, zgrid);
        // printf("x_new y_new z_new: %d %d %d\n", xgrid_new, ygrid_new, zgrid_new);

        double phiA_new, phiB_new;
        double phiA_old_from, phiB_old_from, phiA_old_to, phiB_old_to;
        double phiA_new_from, phiB_new_from, phiA_new_to, phiB_new_to;
        double enrg_local_old, enrg_local_new;
        double enrg_bond_old, enrg_bond_new;
        double enrg_diff;
        int iaccept;

        // Update the new density, ONLY WHEN THE GRID CHANGES!
        if (xgrid_new != xgrid || ygrid_new != xgrid || zgrid_new != zgrid)
        {
            // record the old density of the grid that this particle is:
            // moving from (_from) and moving to (_to)
            phiA_old_from = density_grids[0][xgrid][ygrid][zgrid];
            phiB_old_from = density_grids[1][xgrid][ygrid][zgrid];
            phiA_old_to = density_grids[0][xgrid_new][ygrid_new][zgrid_new];
            phiB_old_to = density_grids[1][xgrid_new][ygrid_new][zgrid_new];

            enrg_local_old = calc_singlemesh_energy(phiA_old_from, phiB_old_from, chi, kappa, density) +
                             calc_singlemesh_energy(phiA_old_to, phiB_old_to, chi, kappa, density);
            if (atom_type == 0)
            {
                phiA_new_from = phiA_old_from - 1 / rho_norm;
                phiA_new_to = phiA_old_to + 1 / rho_norm;
                phiB_new_from = phiB_old_from;
                phiB_new_to = phiB_old_to;
            }
            else
            {
                phiB_new_from = phiB_old_from - 1 / rho_norm;
                phiB_new_to = phiB_old_to + 1 / rho_norm;
                phiA_new_from = phiA_old_from;
                phiA_new_to = phiA_old_to;
            }
            // calculate the single mesh energy
            enrg_local_new = calc_singlemesh_energy(phiA_new_from, phiB_new_from, chi, kappa, density) +
                             calc_singlemesh_energy(phiA_new_to, phiB_new_to, chi, kappa, density);
            // Calculate bond energy
            enrg_bond_old = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);

            // move the particle to the new position and calculate the new bond energy
            coords[iatom][0] = rx_new;
            coords[iatom][1] = ry_new;
            coords[iatom][2] = rz_new;
            enrg_bond_new = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);

            enrg_diff = (enrg_bond_new - enrg_bond_old);

            // Use metropolis criteria to determine whether to accept the move
            iaccept = metro_crit(enrg_diff, idum);
            // For debug:
            // printf("phiA_old_from = %lf, phiA_old_to = %lf\n", phiA_old_from, phiA_old_to);
            // printf("phiB_old_from = %lf, phiB_old_to = %lf\n", phiB_old_from, phiB_old_to);
            // printf("enrg_local_old = %lf, enrg_local_new = %lf\n, \t enrg_bond_old = %lf, enrg_bond_new = %lf\n, \t enrg_diff = %lf, accept = %d\n", enrg_local_old, enrg_local_new, enrg_bond_old, enrg_bond_new, enrg_diff, iaccept);
            *naccept += iaccept;
            if (iaccept == 0)
            {
                // move back the particle if not accepted
                coords[iatom][0] = rx;
                coords[iatom][1] = ry;
                coords[iatom][2] = rz;
            }
            else
            {
                // update the density profile if accepted
                if (atom_type == 0)
                {
                    density_grids[0][xgrid][ygrid][zgrid] = phiA_new_from;
                    density_grids[0][xgrid_new][ygrid_new][zgrid_new] = phiA_new_to;
                }
                else
                {
                    density_grids[1][xgrid][ygrid][zgrid] = phiB_new_from;
                    density_grids[1][xgrid_new][ygrid_new][zgrid_new] = phiB_new_to;
                }
            }
        }
        // only calculate bond energy if the grid is not changed
        else
        {
            // Calculate bond energy
            enrg_bond_old = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);

            // move the particle to the new position and calculate the new bond energy
            coords[iatom][0] = rx_new;
            coords[iatom][1] = ry_new;
            coords[iatom][2] = rz_new;
            enrg_bond_new = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);
            enrg_diff = enrg_bond_new - enrg_bond_old;

            // Use metropolis criteria to determine whether to accept the move
            iaccept = metro_crit(enrg_diff, idum);
            *naccept += iaccept;
            if (iaccept == 0)
            {
                // move back the particle if not accepted
                coords[iatom][0] = rx;
                coords[iatom][1] = ry;
                coords[iatom][2] = rz;
            }
        }
    }

    if (pm_type == 1)
    {
        double rx, ry, rz;
        int xgrid, ygrid, zgrid;
        int atom_type;
        double kappaN, pre_factor;        
        double xc, yc, zc;
        double dx, dy, dz;
        double tx, ty, tz;
        double dx2, dy2, dz2;
        double tx2, ty2, tz2;
        int xneigh, yneigh, zneigh;

        atom_type = particle_type[iatom];
        rx = coords[iatom][0];
        ry = coords[iatom][1];
        rz = coords[iatom][2];

        // identify the grid that the particle is in and the corresponding grid center;
        // taking a grid shift into account
        xgrid = rx / grid_size + grid_shift[0];
        xc = ((double)xgrid + 0.5 - grid_shift[0]) * grid_size;
        xgrid = mod(xgrid, maxsite_1d);
        ygrid = ry / grid_size + grid_shift[1];
        yc = ((double)ygrid + 0.5 - grid_shift[1]) * grid_size;
        ygrid = mod(ygrid, maxsite_1d);
        zgrid = rz / grid_size + grid_shift[2];
        zc = ((double)zgrid + 0.5 - grid_shift[2]) * grid_size;
        zgrid = mod(zgrid, maxsite_1d);

        dx = (rx - xc) / grid_size;
        dy = (ry - yc) / grid_size;
        dz = (rz - zc) / grid_size;

        tx = 1 - fabs(dx);
        ty = 1 - fabs(dy);
        tz = 1 - fabs(dz);

        xneigh = xgrid + ((dx > 0) ? 1 : -1);
        xneigh = mod(xneigh, maxsite_1d);
        yneigh = ygrid + ((dy > 0) ? 1 : -1);
        yneigh = mod(yneigh, maxsite_1d);
        zneigh = zgrid + ((dz > 0) ? 1 : -1);
        zneigh = mod(zneigh, maxsite_1d);
        //if(xgrid < 0 || ygrid < 0 || zgrid < 0 || xneigh < 0 || yneigh < 0 || zneigh < 0)
        //{
        //    printf("xgrid = %d, ygrid = %d, zgrid = %d\n", xgrid, ygrid, zgrid);
        //    printf("xneigh = %d, yneigh = %d, zneigh = %d\n", xneigh, yneigh, zneigh);
        //}
        //if(dx <0 || dy < 0 || dz < 0 || tx < 0 || ty < 0 || tz < 0)
        //{
        //    printf("dx = %f, dy = %f, dz = %f\n", dx, dy, dz);
        //    printf("tx = %f, ty = %f, tz = %f\n", tx, ty, tz);
        //}
        // Randomly displace the particle
        double theta, phi;
        double rx_new, ry_new, rz_new;
        int xgrid_new, ygrid_new, zgrid_new, xneigh_new, yneigh_new, zneigh_new;

        theta = ran1(idum) * 2.0 * M_PI;
        // randomly get cos() value ranging from -1 to 1 to set 2Pi rotation angle from z axis
        phi = acos(2.0 * ran1(idum) - 1.0);
        // move the particle by a vector with length displ
        rx_new = rx + displ * cos(theta) * sin(phi);
        ry_new = ry + displ * sin(theta) * sin(phi);
        rz_new = rz + displ * cos(phi);
        // wrap around the periodic boundary condition
        rx_new = periodic(rx_new, box_size);
        ry_new = periodic(ry_new, box_size);
        rz_new = periodic(rz_new, box_size);

        // identify the grid that the particle is in and the corresponding grid center;
        xgrid_new = rx_new / grid_size + grid_shift[0];
        xc = ((double)xgrid_new + 0.5 - grid_shift[0]) * grid_size;
        xgrid_new = mod(xgrid_new, maxsite_1d);
        ygrid_new = ry_new / grid_size + grid_shift[1];
        yc = ((double)ygrid_new + 0.5 - grid_shift[1]) * grid_size;
        ygrid_new = mod(ygrid_new, maxsite_1d);
        zgrid_new = rz_new / grid_size + grid_shift[2];
        zc = ((double)zgrid_new + 0.5 - grid_shift[2]) * grid_size;
        zgrid_new = mod(zgrid_new, maxsite_1d);

        dx2 = (rx_new - xc) / grid_size;
        dy2 = (ry_new - yc) / grid_size;
        dz2 = (rz_new - zc) / grid_size;

        tx2 = 1 - fabs(dx2);
        ty2 = 1 - fabs(dy2);
        tz2 = 1 - fabs(dz2);

        xneigh_new = xgrid_new + ((dx2 > 0) ? 1 : -1);
        xneigh_new = mod(xneigh_new, maxsite_1d);
        yneigh_new = ygrid_new + ((dy2 > 0) ? 1 : -1);
        yneigh_new = mod(yneigh_new, maxsite_1d);
        zneigh_new = zgrid_new + ((dz2 > 0) ? 1 : -1);
        zneigh_new = mod(zneigh_new, maxsite_1d);

        double phiA_new, phiB_new;
        double phiA_old_from, phiB_old_from, phiA_old_to, phiB_old_to;
        double phiA_new_from, phiB_new_from, phiA_new_to, phiB_new_to;
        double enrg_local_old, enrg_local_new;
        double enrg_bond_old, enrg_bond_new;
        double enrg_diff;
        int iaccept;
        double phiA_old_from1, phiA_old_from2, phiA_old_from3, phiA_old_from4;
        double phiA_old_from5, phiA_old_from6, phiA_old_from7, phiA_old_from8;
        double phiA_old_to1, phiA_old_to2, phiA_old_to3, phiA_old_to4;
        double phiA_old_to5, phiA_old_to6, phiA_old_to7, phiA_old_to8;
        double phiB_old_from1, phiB_old_from2, phiB_old_from3, phiB_old_from4;
        double phiB_old_from5, phiB_old_from6, phiB_old_from7, phiB_old_from8;
        double phiB_old_to1, phiB_old_to2, phiB_old_to3, phiB_old_to4;
        double phiB_old_to5, phiB_old_to6, phiB_old_to7, phiB_old_to8;
        double phiA_new_from1, phiA_new_from2, phiA_new_from3, phiA_new_from4;
        double phiA_new_from5, phiA_new_from6, phiA_new_from7, phiA_new_from8;
        double phiA_new_to1, phiA_new_to2, phiA_new_to3, phiA_new_to4;
        double phiA_new_to5, phiA_new_to6, phiA_new_to7, phiA_new_to8;
        double phiB_new_from1, phiB_new_from2, phiB_new_from3, phiB_new_from4;
        double phiB_new_from5, phiB_new_from6, phiB_new_from7, phiB_new_from8;
        double phiB_new_to1, phiB_new_to2, phiB_new_to3, phiB_new_to4;
        double phiB_new_to5, phiB_new_to6, phiB_new_to7, phiB_new_to8;

        //Update new density when the grid changes
        if(xgrid_new == xgrid && ygrid_new == ygrid && zgrid_new == zgrid && xneigh == xneigh_new && yneigh == yneigh_new && zneigh == zneigh_new)
        {
            //calculate bond energy
            enrg_bond_old = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);
            coords[iatom][0] = rx_new;
            coords[iatom][1] = ry_new;
            coords[iatom][2] = rz_new;
            enrg_bond_new = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);
            //calculate the energy difference
            enrg_diff = enrg_bond_new - enrg_bond_old;        
            // Use metroplis algorithm to decide whether to accept the move
            iaccept = metro_crit(enrg_diff, idum);
            *naccept += iaccept;
            if(iaccept == 0)
            {
                coords[iatom][0] = rx;
                coords[iatom][1] = ry;
                coords[iatom][2] = rz;
            }
        }
        else
        {
            phiA_old_from1 = density_grids[0][xgrid][ygrid][zgrid];
            phiA_old_from2 = density_grids[0][xneigh][ygrid][zgrid];
            phiA_old_from3 = density_grids[0][xgrid][yneigh][zgrid];
            phiA_old_from4 = density_grids[0][xgrid][ygrid][zneigh];
            phiA_old_from5 = density_grids[0][xneigh][yneigh][zgrid];
            phiA_old_from6 = density_grids[0][xneigh][ygrid][zneigh];
            phiA_old_from7 = density_grids[0][xgrid][yneigh][zneigh];
            phiA_old_from8 = density_grids[0][xneigh][yneigh][zneigh];

            phiA_old_to1 = density_grids[0][xgrid_new][ygrid_new][zgrid_new];
            phiA_old_to2 = density_grids[0][xneigh_new][ygrid_new][zgrid_new];
            phiA_old_to3 = density_grids[0][xgrid_new][yneigh_new][zgrid_new];
            phiA_old_to4 = density_grids[0][xgrid_new][ygrid_new][zneigh_new];
            phiA_old_to5 = density_grids[0][xneigh_new][yneigh_new][zgrid_new];
            phiA_old_to6 = density_grids[0][xneigh_new][ygrid_new][zneigh_new];
            phiA_old_to7 = density_grids[0][xgrid_new][yneigh_new][zneigh_new];
            phiA_old_to8 = density_grids[0][xneigh_new][yneigh_new][zneigh_new];

            phiB_old_from1 = density_grids[1][xgrid][ygrid][zgrid];
            phiB_old_from2 = density_grids[1][xneigh][ygrid][zgrid];
            phiB_old_from3 = density_grids[1][xgrid][yneigh][zgrid];
            phiB_old_from4 = density_grids[1][xgrid][ygrid][zneigh];
            phiB_old_from5 = density_grids[1][xneigh][yneigh][zgrid];
            phiB_old_from6 = density_grids[1][xneigh][ygrid][zneigh];
            phiB_old_from7 = density_grids[1][xgrid][yneigh][zneigh];
            phiB_old_from8 = density_grids[1][xneigh][yneigh][zneigh];

            phiB_old_to1 = density_grids[1][xgrid_new][ygrid_new][zgrid_new];
            phiB_old_to2 = density_grids[1][xneigh_new][ygrid_new][zgrid_new];
            phiB_old_to3 = density_grids[1][xgrid_new][yneigh_new][zgrid_new];
            phiB_old_to4 = density_grids[1][xgrid_new][ygrid_new][zneigh_new];
            phiB_old_to5 = density_grids[1][xneigh_new][yneigh_new][zgrid_new];
            phiB_old_to6 = density_grids[1][xneigh_new][ygrid_new][zneigh_new];
            phiB_old_to7 = density_grids[1][xgrid_new][yneigh_new][zneigh_new];
            phiB_old_to8 = density_grids[1][xneigh_new][yneigh_new][zneigh_new];

            enrg_local_old = calc_singlemesh_energy(phiA_old_from1, phiB_old_from1, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from2, phiB_old_from2, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from3, phiB_old_from3, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from4, phiB_old_from4, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from5, phiB_old_from5, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from6, phiB_old_from6, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from7, phiB_old_from7, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_from8, phiB_old_from8, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to1, phiB_old_to1, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to2, phiB_old_to2, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to3, phiB_old_to3, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to4, phiB_old_to4, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to5, phiB_old_to5, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to6, phiB_old_to6, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to7, phiB_old_to7, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_old_to8, phiB_old_to8, chi, kappa, density);

            //Here must have bug...
            if (atom_type == 0)
            {
                phiA_new_from1 = phiA_old_from1 - tx*ty*tz/rho_norm;
                phiA_new_from2 = phiA_old_from2 - (1-tx)*ty*tz/rho_norm;
                phiA_new_from3 = phiA_old_from3 - tx*(1-ty)*tz/rho_norm;
                phiA_new_from4 = phiA_old_from4 - tx*ty*(1-tz)/rho_norm;
                phiA_new_from5 = phiA_old_from5 - (1-tx)*(1-ty)*tz/rho_norm;
                phiA_new_from6 = phiA_old_from6 - (1-tx)*ty*(1-tz)/rho_norm;
                phiA_new_from7 = phiA_old_from7 - tx*(1-ty)*(1-tz)/rho_norm;
                phiA_new_from8 = phiA_old_from8 - (1-tx)*(1-ty)*(1-tz)/rho_norm;
                
                phiA_new_to1 = phiA_old_to1 + tx*ty*tz/rho_norm;
                phiA_new_to2 = phiA_old_to2 + (1-tx)*ty*tz/rho_norm;
                phiA_new_to3 = phiA_old_to3 + tx*(1-ty)*tz/rho_norm;
                phiA_new_to4 = phiA_old_to4 + tx*ty*(1-tz)/rho_norm;
                phiA_new_to5 = phiA_old_to5 + (1-tx)*(1-ty)*tz/rho_norm;
                phiA_new_to6 = phiA_old_to6 + (1-tx)*ty*(1-tz)/rho_norm;
                phiA_new_to7 = phiA_old_to7 + tx*(1-ty)*(1-tz)/rho_norm;
                phiA_new_to8 = phiA_old_to8 + (1-tx)*(1-ty)*(1-tz)/rho_norm;

                phiB_new_from1 = phiB_old_from1;
                phiB_new_from2 = phiB_old_from2;
                phiB_new_from3 = phiB_old_from3;
                phiB_new_from4 = phiB_old_from4;
                phiB_new_from5 = phiB_old_from5;
                phiB_new_from6 = phiB_old_from6;
                phiB_new_from7 = phiB_old_from7;
                phiB_new_from8 = phiB_old_from8;

                phiB_new_to1 = phiB_old_to1;
                phiB_new_to2 = phiB_old_to2;
                phiB_new_to3 = phiB_old_to3;
                phiB_new_to4 = phiB_old_to4;
                phiB_new_to5 = phiB_old_to5;
                phiB_new_to6 = phiB_old_to6;
                phiB_new_to7 = phiB_old_to7;
                phiB_new_to8 = phiB_old_to8;
            }
            else
            {
                phiB_new_from1 = phiB_old_from1 - tx*ty*tz/rho_norm;
                phiB_new_from2 = phiB_old_from2 - (1-tx)*ty*tz/rho_norm;
                phiB_new_from3 = phiB_old_from3 - tx*(1-ty)*tz/rho_norm;
                phiB_new_from4 = phiB_old_from4 - tx*ty*(1-tz)/rho_norm;
                phiB_new_from5 = phiB_old_from5 - (1-tx)*(1-ty)*tz/rho_norm;
                phiB_new_from6 = phiB_old_from6 - (1-tx)*ty*(1-tz)/rho_norm;
                phiB_new_from7 = phiB_old_from7 - tx*(1-ty)*(1-tz)/rho_norm;
                phiB_new_from8 = phiB_old_from8 - (1-tx)*(1-ty)*(1-tz)/rho_norm;

                phiB_new_to1 = phiB_old_to1 + tx*ty*tz/rho_norm;
                phiB_new_to2 = phiB_old_to2 + (1-tx)*ty*tz/rho_norm;
                phiB_new_to3 = phiB_old_to3 + tx*(1-ty)*tz/rho_norm;
                phiB_new_to4 = phiB_old_to4 + tx*ty*(1-tz)/rho_norm;
                phiB_new_to5 = phiB_old_to5 + (1-tx)*(1-ty)*tz/rho_norm;
                phiB_new_to6 = phiB_old_to6 + (1-tx)*ty*(1-tz)/rho_norm;
                phiB_new_to7 = phiB_old_to7 + tx*(1-ty)*(1-tz)/rho_norm;
                phiB_new_to8 = phiB_old_to8 + (1-tx)*(1-ty)*(1-tz)/rho_norm;

                phiA_new_from1 = phiA_old_from1;
                phiA_new_from2 = phiA_old_from2;
                phiA_new_from3 = phiA_old_from3;
                phiA_new_from4 = phiA_old_from4;
                phiA_new_from5 = phiA_old_from5;
                phiA_new_from6 = phiA_old_from6;
                phiA_new_from7 = phiA_old_from7;
                phiA_new_from8 = phiA_old_from8;

                phiA_new_to1 = phiA_old_to1;
                phiA_new_to2 = phiA_old_to2;
                phiA_new_to3 = phiA_old_to3;
                phiA_new_to4 = phiA_old_to4;
                phiA_new_to5 = phiA_old_to5;
                phiA_new_to6 = phiA_old_to6;
                phiA_new_to7 = phiA_old_to7;
                phiA_new_to8 = phiA_old_to8;
            }
            // calculate the single mesh energy
            enrg_local_new = calc_singlemesh_energy(phiA_new_from1, phiB_new_from1, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from2, phiB_new_from2, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from3, phiB_new_from3, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from4, phiB_new_from4, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from5, phiB_new_from5, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from6, phiB_new_from6, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from7, phiB_new_from7, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_from8, phiB_new_from8, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to1, phiB_new_to1, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to2, phiB_new_to2, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to3, phiB_new_to3, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to4, phiB_new_to4, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to5, phiB_new_to5, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to6, phiB_new_to6, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to7, phiB_new_to7, chi, kappa, density)+
                                calc_singlemesh_energy(phiA_new_to8, phiB_new_to8, chi, kappa, density);
            //printf("%f %f %f %f %f %f %f %f\n",calc_singlemesh_energy(phiA_new_from1, phiB_new_from1, chi, kappa, density), calc_singlemesh_energy(phiA_new_from2, phiB_new_from2, chi, kappa, density), calc_singlemesh_energy(phiA_new_from3, phiB_new_from3, chi, kappa, density), calc_singlemesh_energy(phiA_new_from4, phiB_new_from4, chi, kappa, density), calc_singlemesh_energy(phiA_new_from5, phiB_new_from5, chi, kappa, density), calc_singlemesh_energy(phiA_new_from6, phiB_new_from6, chi, kappa, density), calc_singlemesh_energy(phiA_new_from7, phiB_new_from7, chi, kappa, density), calc_singlemesh_energy(phiA_new_from8, phiB_new_from8, chi, kappa, density));
            //calculate bond energy
            enrg_bond_old = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);
            coords[iatom][0] = rx_new;
            coords[iatom][1] = ry_new;
            coords[iatom][2] = rz_new;
            enrg_bond_new = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);
            //calculate the energy difference
            enrg_diff = (enrg_bond_new - enrg_bond_old);  
            //printf("enrg_diff = %f\n", enrg_diff);      
            // Use metroplis algorithm to decide whether to accept the move
            iaccept = metro_crit(enrg_diff, idum);
            //printf("iaccept = %d\n", iaccept);
            *naccept += iaccept;
            if(iaccept == 0)
            {
                coords[iatom][0] = rx;
                coords[iatom][1] = ry;
                coords[iatom][2] = rz;
            }
            else
            {
                coords[iatom][0] = rx_new;
                coords[iatom][1] = ry_new;
                coords[iatom][2] = rz_new;
                
                if(atom_type == 0)
                {
                    density_grids[0][xgrid][ygrid][zgrid] = phiA_new_from1;
                    density_grids[0][xneigh][ygrid][zgrid] = phiA_new_from2;
                    density_grids[0][xgrid][yneigh][zgrid] = phiA_new_from3;
                    density_grids[0][xgrid][ygrid][zneigh] = phiA_new_from4;
                    density_grids[0][xneigh][yneigh][zgrid] = phiA_new_from5;
                    density_grids[0][xneigh][ygrid][zneigh] = phiA_new_from6;
                    density_grids[0][xgrid][yneigh][zneigh] = phiA_new_from7;
                    density_grids[0][xneigh][yneigh][zneigh] = phiA_new_from8;

                    density_grids[0][xgrid_new][ygrid_new][zgrid_new] = phiA_new_to1;
                    density_grids[0][xneigh_new][ygrid_new][zgrid_new] = phiA_new_to2;
                    density_grids[0][xgrid_new][yneigh_new][zgrid_new] = phiA_new_to3;
                    density_grids[0][xgrid_new][ygrid_new][zneigh_new] = phiA_new_to4;
                    density_grids[0][xneigh_new][yneigh_new][zgrid_new] = phiA_new_to5;
                    density_grids[0][xneigh_new][ygrid_new][zneigh_new] = phiA_new_to6;
                    density_grids[0][xgrid_new][yneigh_new][zneigh_new] = phiA_new_to7;
                    density_grids[0][xneigh_new][yneigh_new][zneigh_new] = phiA_new_to8;
                }
                else
                {
                    density_grids[1][xgrid][ygrid][zgrid] = phiB_new_from1;
                    density_grids[1][xneigh][ygrid][zgrid] = phiB_new_from2;
                    density_grids[1][xgrid][yneigh][zgrid] = phiB_new_from3;
                    density_grids[1][xgrid][ygrid][zneigh] = phiB_new_from4;
                    density_grids[1][xneigh][yneigh][zgrid] = phiB_new_from5;
                    density_grids[1][xneigh][ygrid][zneigh] = phiB_new_from6;
                    density_grids[1][xgrid][yneigh][zneigh] = phiB_new_from7;
                    density_grids[1][xneigh][yneigh][zneigh] = phiB_new_from8;

                    density_grids[1][xgrid_new][ygrid_new][zgrid_new] = phiB_new_to1;
                    density_grids[1][xneigh_new][ygrid_new][zgrid_new] = phiB_new_to2;
                    density_grids[1][xgrid_new][yneigh_new][zgrid_new] = phiB_new_to3;
                    density_grids[1][xgrid_new][ygrid_new][zneigh_new] = phiB_new_to4;
                    density_grids[1][xneigh_new][yneigh_new][zgrid_new] = phiB_new_to5;
                    density_grids[1][xneigh_new][ygrid_new][zneigh_new] = phiB_new_to6;
                    density_grids[1][xgrid_new][yneigh_new][zneigh_new] = phiB_new_to7;
                    density_grids[1][xneigh_new][yneigh_new][zneigh_new] = phiB_new_to8;

                }
            
            }
        }

    }
}
