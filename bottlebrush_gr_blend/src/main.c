#include "header.h"

int main(int argc, const char *argv[])
{
    long long int natom_perchain1, natom_perchain2, natom_total1, natom_total2, natom1, natom2, Nbb1, Nbb2, Nsc1, Nsc2, NISC1, NISC2, bbperISC1, bbperISC2, f_branch1, f_branch2, natom_total;
    int nchain_total, nbonds_perchain1, nbonds_perchain2, nbonds_total, nint, nchain1, nchain2;

    int pm_type, cfg, dint;
    double kappa, kspring, chi, bondlen, bondlen_sqr, Nbar;
    double f_A;
    double grid_size;
    int maxsite_1d;

    int save_every, tot_mc_cycles, equil_cycles, grid_shift_every;
    double ave_bondlen, reedsq_bb, reedsq_sc, rgsq_bb, rgsq_sc;
    double reedsq_bb1, reedsq_bb2, reedsq_sc1, reedsq_sc2, rgsq_bb1, rgsq_bb2, rgsq_sc1, rgsq_sc2;

    time_t start_t, current_t;
    double diff_t;
    int diff_days, diff_hours, diff_minutes, diff_seconds;
    time(&start_t);

    // File Input
    FILE *inputfile;
    inputfile = fopen("Input.txt", "r");
    fscanf(inputfile, "Nbb1 = %lld\n", &Nbb1);           // # of backbone atoms per chain
    fscanf(inputfile, "Nsc1 = %lld\n", &Nsc1);           // # of backbone atoms per chain
    fscanf(inputfile, "f_branch1 = %lld\n", &f_branch1); // # of backbone atoms per chain
    fscanf(inputfile, "Nbb2 = %lld\n", &Nbb2);           // # of backbone atoms per chain
    fscanf(inputfile, "Nsc2 = %lld\n", &Nsc2);           // # of backbone atoms per chain
    fscanf(inputfile, "f_branch2 = %lld\n", &f_branch2); // # of backbone atoms per chain
    fscanf(inputfile, "frac_A = %lf\n", &f_A);            // volume fraction of A
    fscanf(inputfile, "pm_type = %d\n", &pm_type);  // PM type (0- or 1-order)
    fscanf(inputfile, "chi = %lf\n", &chi);         // F-H parameter
    fscanf(inputfile, "kappa = %lf\n", &kappa);     // imcompressible parameter
    fscanf(inputfile, "kspring = %lf\n", &kspring); // spring constant
    fscanf(inputfile, "grid_size(Re) = %lf\n", &grid_size);          // mesh size
    fscanf(inputfile, "maxsite_1d = %d\n", &maxsite_1d);             // number of grids per direction
    fscanf(inputfile, "save_every = %d\n", &save_every);             // saving interval
    fscanf(inputfile, "tot_mc_cycles = %d\n", &tot_mc_cycles);       // # of steps (exclude equil. steps)
    fscanf(inputfile, "equil_cycles = %d\n", &equil_cycles);         // # of equil. steps
    fscanf(inputfile, "grid_shift_every = %d\n", &grid_shift_every); // frequency of shift grid
    fscanf(inputfile, "NISC1 = %d\n", &NISC1);
    fscanf(inputfile, "NISC2 = %d\n", &NISC2);
    fscanf(inputfile, "bbperISC1 = %d\n", &bbperISC1);
    fscanf(inputfile, "bbperISC2 = %d\n", &bbperISC2);  
    fscanf(inputfile, "cfg = %d\n", &cfg);                 

    fclose(inputfile);

    double density, rho_norm;
    double box_size, box_volume;

    int nsites;
    bondlen=1.0;
    bondlen_sqr = bondlen * bondlen;
    natom_perchain1 = Nbb1 * (f_branch1 * Nsc1 + 1);
    natom_perchain2 = Nbb2 * (f_branch2 * Nsc2 + 1);
    double Re, ReReRe;
    Re = sqrt((20 - 1) * bondlen_sqr);
    ReReRe = Re * Re * Re;
    grid_size = grid_size * Re;
    density=5*20/ReReRe;
    nint = density * grid_size * grid_size * grid_size;
    // IMPORTANT NOTE:
    // "density" refers to the number density of bead, whereas
    // "rho_norm" refers to the normalized density in perfectly homogeneous melt
    // Here, assume nint = 30, density = nint / deltaL**3
    //density = 30 / (grid_size * grid_size * grid_size);

    // spring constant is 1.5 / b^2
    kspring = kspring / bondlen_sqr;
    box_size = maxsite_1d * grid_size;
    box_volume = box_size * box_size * box_size;
    nsites = maxsite_1d * maxsite_1d * maxsite_1d;
    nchain1 = density * box_volume*f_A / natom_perchain1;
    nchain2 = density * box_volume*(1-f_A) / natom_perchain2;
    nchain_total = nchain1 + nchain2;  
    natom1 = nchain1 * natom_perchain1;   
    natom2 = nchain2 * natom_perchain2;   
    natom_total = natom1+natom2; 
    nbonds_perchain1 = natom_perchain1 - 1;  
    nbonds_perchain2 = natom_perchain2 - 1;  
    nbonds_total = natom_total - nchain1 - nchain2;   
    rho_norm = (double)natom_total / nsites; // normalized density homogeneous melt blend

    FILE *outputfile;
    outputfile = fopen("output.txt", "w");
    fprintf(outputfile, "nchain1 = %d\n", nchain1);
    fprintf(outputfile, "nchain2 = %d\n", nchain2);
    fprintf(outputfile, "natom_perchain1 = %lld\n", natom_perchain1);
    fprintf(outputfile, "natom_perchain2 = %lld\n", natom_perchain2);
    fprintf(outputfile, "natom_total = %lld\n", natom_total);
    fprintf(outputfile, "nbonds_total = %d\n", nbonds_total);
    fprintf(outputfile, "nint = %d\n", nint);
    fprintf(outputfile, "density = %lf\n", density);
    fprintf(outputfile, "rho_norm = %lf\n", rho_norm);
    fprintf(outputfile, "bondlen = %lf\n", bondlen);
    fprintf(outputfile, "chi = %lf\n", chi);
    fprintf(outputfile, "kappa = %lf\n", kappa);
    fprintf(outputfile, "kspring = %lf\n", kspring);
    fprintf(outputfile, "box_size = %lf\n", box_size);
    fprintf(outputfile, "grid_size = %lf\n", grid_size);
    fprintf(outputfile, "maxsite_1d = %d\n", maxsite_1d);
    fprintf(outputfile, "nsites = %d\n", nsites);

    fclose(outputfile);

    // Memory allocation
    // type of particle, 1d array (natoms)
    int *particle_type = calloc(natom_total, sizeof(int));
    // Type of bond 
    int *bond_type = calloc(nbonds_total, sizeof(int));
    // molecule ID, 1d array (natoms)
    int *molecule_id = calloc(natom_total, sizeof(int));
    // coordinates of particles, 2d array (natoms, 3)
    double **coords;
    coords = (double **)calloc(natom_total, sizeof(double *));
    for (int i = 0; i < natom_total; ++i)
    {
        coords[i] = (double *)calloc(3, sizeof(double));
    }
    coords[0][0] = 0;
    //bond table, 2d array (nbonds, 3)
    int **bonds;
    bonds = (int **)calloc(nbonds_total, sizeof(int *));
    for (int i = 0; i < nbonds_total; ++i)
    {
        bonds[i] = (int *)calloc(2, sizeof(int));
    }
    //atom bond list, each entry lists the atoms that the given atom is bonded with
    int **atom_bond_list;
    atom_bond_list = (int **)calloc(natom_total, sizeof(int *));
    for (int i = 0; i < natom_total; ++i)
    {
        // Assuming maximum 10 bonds per atom!!!
        atom_bond_list[i] = (int *)calloc(12, sizeof(int));
    }
    // Then initialize this list to -1 to make sure check can be performed
    for (int i = 0; i < natom_total; ++i)
    {
        // Assuming maximum 10 bonds per atom!!!
        for (int j = 0; j < 12; ++j)
        {
            atom_bond_list[i][j] = -1;
        }
    }
    // printf("done\n");
   // density profile, 4d array (num_atom_types, num_grid, num_grid, num_grid)
    double ****density_grids = calloc(2, sizeof(double ***));
    for (int itype = 0; itype < 2; ++itype)
    {
        density_grids[itype] = calloc(maxsite_1d, sizeof(double **));
        for (int i = 0; i < maxsite_1d; ++i)
        {
            density_grids[itype][i] = calloc(maxsite_1d, sizeof(double *));
            for (int j = 0; j < maxsite_1d; ++j)
            {
                density_grids[itype][i][j] = calloc(maxsite_1d, sizeof(double));
            }
        }
    }
    double *avgRDF_AA = calloc(maxsite_1d*10, sizeof(double));
    double *avgRDF_BB = calloc(maxsite_1d*10, sizeof(double));
    double *avgRDF_AB = calloc(maxsite_1d*10, sizeof(double));
    double *grAA = calloc(maxsite_1d*10, sizeof(double));
    double *grBB = calloc(maxsite_1d*10, sizeof(double));
    double *grAB = calloc(maxsite_1d*10, sizeof(double)); 
    // printf("done\n");

    // Initialize random number generator
    long *idum = malloc(sizeof(long));
    *idum = initRan();

    // Initialize polymers
    int start_idx = 0;
    int mol_start_idx = 0;
    if(cfg==0)
    {
        initialize_brushblend(nchain1, nchain2, 0, 1, 1, 1, Nbb1, Nbb2, Nsc1, Nsc2, f_branch1, f_branch2, box_size, start_idx, mol_start_idx, particle_type, bond_type, molecule_id,
                                coords, bonds, atom_bond_list, idum);
        // Write the initial configuration to a file
        FILE *init_config;
        init_config = fopen("ini.cfg", "w");

        fprintf(init_config, "#Initial configuration\n");
        fprintf(init_config, "%lld atoms\n", natom_total);
        fprintf(init_config, "2 atom types\n");
        fprintf(init_config, "%d bonds\n", nbonds_total);
        fprintf(init_config, "2 bond types\n");
        fprintf(init_config, "0.0 %lf xlo xhi\n", box_size);
        fprintf(init_config, "0.0 %lf ylo yhi\n", box_size);
        fprintf(init_config, "0.0 %lf zlo zhi\n", box_size);
        fprintf(init_config, "0 0 0 xy xz yz\n");
        fprintf(init_config, "\n");
        fprintf(init_config, "Atoms # bond\n");
        fprintf(init_config, "\n");

        for (int i = 0; i < natom_total; ++i)
        {
            fprintf(init_config, "%d %d %d %f %f %f\n", i + 1, molecule_id[i], particle_type[i] + 1,
                    coords[i][0], coords[i][1], coords[i][2]);
        }

        fprintf(init_config, "\n");
        fprintf(init_config, "Bonds\n");
        fprintf(init_config, "\n");

        for (int i = 0; i < nbonds_total; ++i)
        {
            fprintf(init_config, "%d %d %d %d\n", i + 1, bond_type[i], bonds[i][0] + 1, bonds[i][1] + 1);
        }
        fclose(init_config);
    }
    else if(cfg == 1)
    {
        initialize_brushblend(nchain1, nchain2, 0, 1, 1, 1, Nbb1, Nbb2, Nsc1, Nsc2, f_branch1, f_branch2, box_size, start_idx, mol_start_idx, particle_type, bond_type, molecule_id,
                                coords, bonds, atom_bond_list, idum);
        // Read the initial configuration from a file
        FILE *init_config;
        init_config = fopen("coord.xyz", "r"); // has only coordinates
        // fscanf(init_config, "#Initial configuration\n");
        // fscanf(init_config, "%lld atoms\n", &natom_total);
        // fscanf(init_config, "2 atom types\n");
        // fscanf(init_config, "%d bonds\n", &nbonds_total);
        // fscanf(init_config, "2 bond types\n");
        // fscanf(init_config, "0.0 %lf xlo xhi\n", &box_size);
        // fscanf(init_config, "0.0 %lf ylo yhi\n", &box_size);
        // fscanf(init_config, "0.0 %lf zlo zhi\n", &box_size);
        // fscanf(init_config, "0 0 0 xy xz yz\n");
        // fscanf(init_config, "\n");
        // fscanf(init_config, "Atoms # bond\n");
        // fscanf(init_config, "\n");

        for (int i = 0; i < natom_total; ++i)
        {
            fscanf(init_config, "%lf %lf %lf\n", &coords[i][0], &coords[i][1], &coords[i][2]);
        }
        // for(int i=0; i<natom_total; i++)
        // {
        //     particle_type[i] -= 1;
        // }

        // fscanf(init_config, "\n");
        // fscanf(init_config, "Bonds\n");
        // fscanf(init_config, "\n");

        // for (int i = 0; i < nbonds_total; ++i)
        // {
        //     fscanf(init_config, "%d %d %d %d\n", &dint, &bond_type[i], &bonds[i][0], &bonds[i][1]);
        // }
        // for(int i=0; i<nbonds_total; i++)
        // {
        //     bonds[i][0] -= 1;
        //     bonds[i][1] -= 1;
        // }
        fclose(init_config);
    }    
    // printf("done\n");

    // Get the initial density profile, with grid shift = 0
    double grid_shift[3];
    grid_shift[0] = 0;
    grid_shift[1] = 0;
    grid_shift[2] = 0;
    get_mesh_density(natom_total, grid_size, maxsite_1d, pm_type, rho_norm,
                     grid_shift, coords, particle_type, density_grids);

    // equilibration run
    // Write the thermo info to a file
    double mesh_energy, bond_energy, total_energy;

    FILE *thermo;
    char *thermo_filename = malloc(sizeof(char) * FILENAME_MAX);
    sprintf(thermo_filename, "thermo_blend.txt");
    thermo = fopen(thermo_filename, "w");
    fprintf(thermo, "t\tmesh_energy\tbond_energy\ttotal_energy\tacceptance_rate\tave_bondlen\treedsq_bb\treed_sc\trgsq_bb\trgsq_sc\n");
    fclose(thermo);
    // // Write the density info to a file
    // FILE *density_profile;
    // char *density_filename = malloc(sizeof(char) * FILENAME_MAX);
    // sprintf(density_filename, "dens_profile_blend.txt");
    // density_profile = fopen(density_filename, "w");
    // fprintf(density_profile, "type\tx\ty\tz\tdensity\n");
    // fclose(density_profile);

    // Wrtie gr info to a file
    FILE *gr_profile;
    char *gr_filename = malloc(sizeof(char) * FILENAME_MAX);
    sprintf(gr_filename, "gr_profile_blend.txt");
    gr_profile = fopen(gr_filename, "w");

    // Write traj
    FILE *traj;
    char *traj_name = malloc(sizeof(char) * FILENAME_MAX);
    sprintf(traj_name, "traj_blend.xyz");
    traj = fopen(traj_name, "w");
    fclose(traj);

    // Dummy coords file
    FILE *dump1, *dump2, *dump3;
    char *dump1_name = malloc(sizeof(char) * FILENAME_MAX);
    char *dump2_name = malloc(sizeof(char) * FILENAME_MAX);
    char *dump3_name = malloc(sizeof(char) * FILENAME_MAX);
    sprintf(dump1_name, "dump1");
    sprintf(dump2_name, "dump2");
    sprintf(dump3_name, "dump3");

    // Write the bond info to dump2
    dump2 = fopen(dump2_name, "w");
    for (int i = 0; i < nbonds_total; ++i)
    {
        fprintf(dump2, "%d %d\n", bonds[i][0], bonds[i][1]);
    }
    fclose(dump2);

    // Write the atom_bond_list info to dump3
    dump3 = fopen(dump3_name, "w");
    for (int i = 0; i < natom_total; ++i)
    {
        for (int j = 0; j < 12; ++j)
        {
            fprintf(dump3, "%d ", atom_bond_list[i][j]);
        }
        fprintf(dump3, "\n");
    }
    fclose(dump3);
    // printf("done\n");

    thermo = fopen(thermo_filename, "a"); 
    traj = fopen(traj_name, "a");
    //density_profile = fopen(density_filename, "a");
    unsigned long long t;
    unsigned long long naccept = 0;
    unsigned long long nattempt = 0;
    double acc_rate = 0.0;
    
    if(cfg==0)
    {
        for (t = 0; t <= equil_cycles; ++t)
        {
            // use a new grid_shift every so often
            if (t % grid_shift_every == 0)
            {
                // reset density grids
                for (int itype = 0; itype < 2; ++itype)
                {
                    for (int i = 0; i < maxsite_1d; ++i)
                    {
                        for (int j = 0; j < maxsite_1d; ++j)
                        {
                            memset(density_grids[itype][i][j], 0, maxsite_1d * sizeof(double));
                        }
                    }
                }
                // move by atmost half the grid size is enough to mimic a random grid discretization
                grid_shift[0] = 0.5 * ran1(idum);
                grid_shift[1] = 0.5 * ran1(idum);
                grid_shift[2] = 0.5 * ran1(idum);
                get_mesh_density(natom_total, grid_size, maxsite_1d, pm_type, rho_norm,
                                grid_shift, coords, particle_type, density_grids);
            }
            // printf("done\n");

            // MC atom displacement move for each particle in the system
            for (int i = 0; i < natom_total; ++i)
            {
                mc_atom_displ_bond(i, bondlen, pm_type, &naccept, density_grids, density,
                            particle_type, coords, atom_bond_list, box_size, grid_size, rho_norm, maxsite_1d, chi, kappa,
                            kspring, natom_total, grid_shift, idum);
                nattempt += 1;
            }

            // save thermo info every so often
            if (t % save_every == 0)
            {
                // printf("natom1: %d\n", natom1);
                acc_rate = (double)naccept / nattempt;
                mesh_energy = calc_total_particlemesh_energy(density_grids, chi, kappa, density, maxsite_1d);
                bond_energy = calc_total_bond_energy_harmonic(coords, bonds, nbonds_total, kspring, box_size);
                total_energy = mesh_energy + bond_energy;
                ave_bondlen = calc_ave_bondlen(coords, bonds, nbonds_total, box_size);
                reedsq_bb1 = calc_reedsq_bb(coords, Nbb1, nchain1, natom_perchain1, box_size, 0);
                reedsq_bb2 = calc_reedsq_bb(coords, Nbb2, nchain2, natom_perchain2, box_size, natom1);
                reedsq_sc1 = calc_reedsq_sc(coords, Nbb1, Nsc1, f_branch1, nchain1, box_size,0);
                reedsq_sc2 = calc_reedsq_sc(coords, Nbb2, Nsc2, f_branch2, nchain2, box_size, natom1);
                rgsq_bb1 = calc_rgsq_bb(coords, Nbb1, nchain1, natom_perchain1, box_size, 0);
                rgsq_bb2 = calc_rgsq_bb(coords, Nbb2, nchain2, natom_perchain2, box_size, natom1);
                rgsq_sc1 = calc_rgsq_sc(coords, Nbb1, Nsc1, f_branch1, nchain1, box_size, 0);
                rgsq_sc2 = calc_rgsq_sc(coords, Nbb2, Nsc2, f_branch2, nchain2, box_size, natom1);
                reedsq_bb = 0.5*(reedsq_bb1 + reedsq_bb2);
                reedsq_sc = 0.5*(reedsq_sc1 + reedsq_sc2);
                rgsq_bb = 0.5*(rgsq_bb1 + rgsq_bb2);
                rgsq_sc = 0.5*(rgsq_sc1 + rgsq_sc2);
                fprintf(thermo, "%llu\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, mesh_energy, bond_energy, total_energy, acc_rate,
                        ave_bondlen, reedsq_bb, reedsq_sc, rgsq_bb, rgsq_sc);
                // printf("%lf %lf %lf\n", reedsq_bb1, reedsq_bb2, reedsq_bb);
            }
        
            // Every save_every*500, print the coordinates and free/reallocate memeory
            if(t % save_every == 0)
            {
                // save coord in dump
                dump1 = fopen(dump1_name, "w");
                for (int i = 0; i < natom_total; ++i)
                {
                    fprintf(dump1, "%lf %lf %lf\n", coords[i][0], coords[i][1], coords[i][2]);
                }
                fclose(dump1);
                // printf("done\n");

                //free memory coords, bonds, atom_bond_list
                for (int i = 0; i < natom_total; ++i)
                {
                    free(coords[i]);
                }
                free(coords);
                // printf("done\n");
                for (int i = 0; i < nbonds_total; ++i)
                {
                    free(bonds[i]);
                }
                free(bonds);
                // printf("done\n");
                for (int i = 0; i < natom_total; ++i)
                {
                    free(atom_bond_list[i]);
                }
                free(atom_bond_list);
                // printf("done\n");

                //reallocate memeory
                coords = (double **)calloc(natom_total, sizeof(double *));
                for (int i = 0; i < natom_total; ++i)
                {
                    coords[i] = (double *)calloc(3, sizeof(double));
                }
                coords[0][0] = 0;
                bonds = (int **)calloc(nbonds_total, sizeof(int *));
                for (int i = 0; i < nbonds_total; ++i)
                {
                    bonds[i] = (int *)calloc(2, sizeof(int));
                }
                atom_bond_list = (int **)calloc(natom_total, sizeof(int *));
                for (int i = 0; i < natom_total; ++i)
                {
                    atom_bond_list[i] = (int *)calloc(12, sizeof(int));
                }

                // read the last configuration of traj file
                dump1 = fopen(dump1_name, "r");
                for (int i = 0; i < natom_total; ++i)
                {
                    fscanf(dump1, "%lf %lf %lf\n", &coords[i][0], &coords[i][1], &coords[i][2]);
                }
                fclose(dump1);
                dump2 = fopen(dump2_name, "r");
                for (int i = 0; i < nbonds_total; ++i)
                {
                    fscanf(dump2, "%d %d\n", &bonds[i][0], &bonds[i][1]);
                }
                fclose(dump2);
                dump3 = fopen(dump3_name, "r");
                for (int i = 0; i < natom_total; ++i)
                {
                    for (int j = 0; j < 12; ++j)
                    {
                        fscanf(dump3, "%d ", &atom_bond_list[i][j]);
                    }
                }
                fclose(dump3);
                // printf("done\n");

            }
        }

        time(&current_t);
        diff_t = difftime(current_t, start_t);
        diff_days = diff_t / 86400;
        diff_hours = (diff_t - diff_days * 86400) / 3600;
        diff_minutes = (diff_t - diff_days * 86400 - diff_hours * 3600) / 60;
        diff_seconds = diff_t - diff_days * 86400 - diff_hours * 3600 - diff_minutes * 60;

        printf("Time elapsed: %dd-%dh-%dm-%ds\n", diff_days, diff_hours, diff_minutes, diff_seconds);
        printf("End of Equilibration\n");
        fprintf(thermo, "\nEnd of Equilibration\n\n");
    }
    else{}
    // printf("done\n");


    // production run
    naccept = 0;
    nattempt = 0;
    for (t = 0; t <= tot_mc_cycles; ++t)
    {
        // use a new grid_shift every so often
        if (t % grid_shift_every == 0)
        {
            for (int itype = 0; itype < 2; ++itype)
            {
                for (int i = 0; i < maxsite_1d; ++i)
                {
                    for (int j = 0; j < maxsite_1d; ++j)
                    {
                        memset(density_grids[itype][i][j], 0, maxsite_1d * sizeof(double));
                    }
                }
            }
            grid_shift[0] = 0.5 * ran1(idum);
            grid_shift[1] = 0.5 * ran1(idum);
            grid_shift[2] = 0.5 * ran1(idum);
            get_mesh_density(natom_total, grid_size, maxsite_1d, pm_type, rho_norm,
                             grid_shift, coords, particle_type, density_grids);
        }

        // MC atom displacement move for each particle in the system
        for (int i = 0; i < natom_total; ++i)
        {
            mc_atom_displ(i, bondlen, pm_type, &naccept, density_grids, density,
                          particle_type, coords, atom_bond_list, box_size, grid_size, rho_norm, maxsite_1d, chi, kappa,
                          kspring, natom_total, grid_shift, idum);
            nattempt += 1;
        }

        // Every save_every*500, print the coordinates and free/reallocate memeory
        // if(t % (save_every*5) == 0)
        // {
        //     // save coord in traj
        //     fprintf(traj, "%lld\n%llu 1.500000\n", natom_total, t);
        //     for (int i = 0; i < natom_total; ++i)
        //     {
        //         fprintf(traj, "%d %lf %lf %lf\n", particle_type[i]+1, coords[i][0], coords[i][1], coords[i][2]);
        //     }
        // }
        // printf("done\n");

        if(t % save_every == 0)
        {
            // save coord in dump
            dump1 = fopen(dump1_name, "w");
            for (int i = 0; i < natom_total; ++i)
            {
                fprintf(dump1, "%lf %lf %lf\n", coords[i][0], coords[i][1], coords[i][2]);
            }
            fclose(dump1);

            //free memory coords, bonds, atom_bond_list
            for (int i = 0; i < natom_total; ++i)
            {
                free(coords[i]);
            }
            free(coords);
            for (int i = 0; i < nbonds_total; ++i)
            {
                free(bonds[i]);
            }
            free(bonds);
            for (int i = 0; i < natom_total; ++i)
            {
                free(atom_bond_list[i]);
            }
            free(atom_bond_list);

            //reallocate memeory
            coords = (double **)calloc(natom_total, sizeof(double *));
            for (int i = 0; i < natom_total; ++i)
            {
                coords[i] = (double *)calloc(3, sizeof(double));
            }
            coords[0][0] = 0;
            bonds = (int **)calloc(nbonds_total, sizeof(int *));
            for (int i = 0; i < nbonds_total; ++i)
            {
                bonds[i] = (int *)calloc(2, sizeof(int));
            }
            atom_bond_list = (int **)calloc(natom_total, sizeof(int *));
            for (int i = 0; i < natom_total; ++i)
            {
                atom_bond_list[i] = (int *)calloc(12, sizeof(int));
            }
            // printf("done\n");

            // read the last configuration of traj file
            dump1 = fopen(dump1_name, "r");
            for (int i = 0; i < natom_total; ++i)
            {
                fscanf(dump1, "%lf %lf %lf\n", &coords[i][0], &coords[i][1], &coords[i][2]);
            }
            fclose(dump1);
            dump2 = fopen(dump2_name, "r");
            for (int i = 0; i < nbonds_total; ++i)
            {
                fscanf(dump2, "%d %d\n", &bonds[i][0], &bonds[i][1]);
            }
            fclose(dump2);
            dump3 = fopen(dump3_name, "r");
            for (int i = 0; i < natom_total; ++i)
            {
                for (int j = 0; j < 12; ++j)
                {
                    fscanf(dump3, "%d ", &atom_bond_list[i][j]);
                }
            }
            fclose(dump3);
        }
        //COM and midCOM 
        // if(t % (save_every*500) == 0)
        // {
        //     // save coord in traj
        //     fprintf(traj, "%lld\n%llu 1.500000\n", natom_total, t);
        //     for (int i = 0; i < natom_total; ++i)
        //     {
        //         if(i/natom_perchain1 < nchain1)
        //         {
        //             fprintf(traj, "%d %lf %lf %lf\n", 1, com[i][0], com[i][1], com[i][2]);
        //             fprintf(traj, "%d %lf %lf %lf\n", 1, midcom[i][0], midcom[i][1], midcom[i][2]);
        //         }
        //         else
        //         {
        //             fprintf(traj, "%d %lf %lf %lf\n", 2, com[i][0], com[i][1], com[i][2]);
        //             fprintf(traj, "%d %lf %lf %lf\n", 2, midcom[i][0], midcom[i][1], midcom[i][2]);
        //         }
        //     }
        // }
        // if (t > (tot_mc_cycles - 1000) && t % save_every == 0)
        // {
        //     // density profile
        //     fprintf(density_profile, "t = %llu\n", t);
        //     for (int itype = 0; itype < 2; ++itype)
        //     {
        //         for (int i = 0; i < maxsite_1d; ++i)
        //         {
        //             for (int j = 0; j < maxsite_1d; ++j)
        //             {
        //                 for (int k = 0; k < maxsite_1d; ++k)
        //                 {
        //                     fprintf(density_profile, "%d\t%d\t%d\t%d\t%f\n", itype, i, j, k, density_grids[itype][i][j][k]);
        //                 }
        //             }
        //         }
        //     }

        // gr profile
        calc_gr_ISC(coords, grAA, grBB, grAB, particle_type, Nbb1, Nbb2, Nsc1, Nsc2, f_branch1, f_branch2, nchain1, nchain2, natom_perchain1, natom_perchain2, box_size, NISC1, NISC2, bbperISC1, bbperISC2, maxsite_1d*10);
        for(int i=0; i<maxsite_1d*10; i++)
        {
            avgRDF_AA[i] += grAA[i];
            avgRDF_BB[i] += grBB[i];
            avgRDF_AB[i] += grAB[i];
        }
        if((t % (save_every*10) == 0) && t!=0)
        {
            fprintf(gr_profile, "r\tgrAA\tgrBB\tgrAB\n");
            fprintf(gr_profile, "#t = %llu\n", t);
            for (int i = 1; i < maxsite_1d*10; ++i)
            {
                fprintf(gr_profile, "%lf\t%lf\t%lf\t%lf\n", (double)box_size/(2*maxsite_1d*10)*i, (double)avgRDF_AA[i]/(save_every*10+1), (double)avgRDF_BB[i]/(save_every*10+1), (double)avgRDF_AB[i]/(save_every*10+1));     
                avgRDF_AA[i] = 0;
                avgRDF_BB[i] = 0;
                avgRDF_AB[i] = 0;
            }
        }
        // printf("done\n");

        // save thermo and density profile every so often
        if (t % save_every == 0)
        {
            time(&current_t);
            diff_t = difftime(current_t, start_t);
            diff_days = diff_t / 86400;
            diff_hours = (diff_t - diff_days * 86400) / 3600;
            diff_minutes = (diff_t - diff_days * 86400 - diff_hours * 3600) / 60;
            diff_seconds = diff_t - diff_days * 86400 - diff_hours * 3600 - diff_minutes * 60;

            printf("Timestep %llu finished. Time elapsed: %dd-%dh-%dm-%ds.\n", t, diff_days, diff_hours, diff_minutes, diff_seconds);
            fflush(stdout);

            // thermo
            acc_rate = (double)naccept / nattempt;
            mesh_energy = calc_total_particlemesh_energy(density_grids, chi, kappa, density, maxsite_1d);
            bond_energy = calc_total_bond_energy_harmonic(coords, bonds, nbonds_total, kspring, box_size);
            total_energy = mesh_energy + bond_energy;
            ave_bondlen = calc_ave_bondlen(coords, bonds, nbonds_total, box_size);
            reedsq_bb1 = calc_reedsq_bb(coords, Nbb1, nchain1, natom_perchain1, box_size,0);
            reedsq_bb2 = calc_reedsq_bb(coords, Nbb2, nchain2, natom_perchain2, box_size,natom1);
            reedsq_sc1 = calc_reedsq_sc(coords, Nbb1, Nsc1, f_branch1, nchain1, box_size,0);
            reedsq_sc2 = calc_reedsq_sc(coords, Nbb2, Nsc2, f_branch2, nchain2, box_size,natom1);
            rgsq_bb1 = calc_rgsq_bb(coords, Nbb1, nchain1, natom_perchain1, box_size, 0);
            rgsq_bb2 = calc_rgsq_bb(coords, Nbb2, nchain2, natom_perchain2, box_size, natom1);
            rgsq_sc1 = calc_rgsq_sc(coords, Nbb1, Nsc1, f_branch1, nchain1, box_size, 0);
            rgsq_sc2 = calc_rgsq_sc(coords, Nbb2, Nsc2, f_branch2, nchain2, box_size, natom1);
            reedsq_bb = 0.5*(reedsq_bb1 + reedsq_bb2);
            reedsq_sc = 0.5*(reedsq_sc1 + reedsq_sc2);
            rgsq_bb = 0.5*(rgsq_bb1 + rgsq_bb2);
            rgsq_sc = 0.5*(rgsq_sc1 + rgsq_sc2);
            fprintf(thermo, "%llu\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, mesh_energy, bond_energy, total_energy, acc_rate, ave_bondlen, reedsq_bb, reedsq_sc, rgsq_bb, rgsq_sc);
            //fflush(thermo);
            //printf("t = %llu, mesh_energy = %lf, bond_energy = %lf, total_energy = %lf, acc_rate = %lf, ave_bondlen = %lf, reedsq_bb = %lf, reedsq_sc = %lf, rgsq_bb = %lf, rgsq_sc = %lf\n", t, mesh_energy, bond_energy, total_energy, acc_rate, ave_bondlen, reedsq_bb, reedsq_sc, rgsq_bb, rgsq_sc);
        }
    }
    fclose(traj);
    fclose(thermo);
    fclose(gr_profile);

    // Write the final configuration to a file
    FILE *final_config;
    final_config = fopen("final.cfg", "w");
    fprintf(final_config, "#Final configuration\n");
    fprintf(final_config, "%lld atoms\n", natom_total);
    fprintf(final_config, "2 atom types\n");
    fprintf(final_config, "%d bonds\n", nbonds_total);
    fprintf(final_config, "2 bond types\n");
    fprintf(final_config, "0.0 %lf xlo xhi\n", box_size);
    fprintf(final_config, "0.0 %lf ylo yhi\n", box_size);
    fprintf(final_config, "0.0 %lf zlo zhi\n", box_size);
    fprintf(final_config, "0 0 0 xy xz yz\n");
    fprintf(final_config, "\n");
    fprintf(final_config, "Atoms # bond\n");
    fprintf(final_config, "\n");
    for (int i = 0; i < natom_total; ++i)
    {
        fprintf(final_config, "%d %d %d %f %f %f\n", i + 1, molecule_id[i], particle_type[i] + 1,
                coords[i][0], coords[i][1], coords[i][2]);
    }
    fprintf(final_config, "\n");
    fprintf(final_config, "Bonds\n");
    fprintf(final_config, "\n");
    for (int i = 0; i < nbonds_total; ++i)
    {
        fprintf(final_config, "%d %d %d %d\n", i + 1, bond_type[i], bonds[i][0] + 1, bonds[i][1] + 1);
    }
    fclose(final_config);
    // printf("done\n");

    // Free memory
    free(particle_type); free(bond_type); free(molecule_id); free(idum); free(thermo_filename); free(gr_filename);
    free(avgRDF_AA); free(avgRDF_AB); free(grAA); free(grAB);
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<maxsite_1d; j++)
        {
            free(density_grids[i][j]);
        }
        free(density_grids[i]);
    }
    free(density_grids);
    for(int i=0; i<natom_total; i++)
    {
        free(coords[i]);
        free(atom_bond_list[i]);
    }
    free(atom_bond_list);free(coords);
    for(int i=0; i<nbonds_total; i++)
    {
        free(bonds[i]);
    }
    free(bonds);
    // printf("done\n");

    //remove dump files
    remove(dump1_name);
    remove(dump2_name);
    remove(dump3_name);
}
