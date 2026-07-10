#include "header.h"
#include "math.h"
// initialization of linear polymer
void initialize_polymer(long long int nchain, long long int natom_perchain, int itype, double bb, double box_size,
                        int start_idx, int mol_start_idx, int *particle_type, int *molecule_id,
                        double **coords, int **bonds, int **atom_bond_list, long *idum)
{
    double theta, phi;
    int index;
    double b = sqrt(bb);
    // double b = 1;
    int ibond;
    ibond = start_idx - mol_start_idx;
    for (int i = 0; i < nchain; ++i)
    {
        // Randomizing positions/locations of the first bead(=index 0, 1N, 2N-th...) of each polymer
        index = i * natom_perchain + start_idx;
        coords[index][0] = ran1(idum) * box_size;
        coords[index][1] = ran1(idum) * box_size;
        coords[index][2] = ran1(idum) * box_size;
        particle_type[index] = itype;
        molecule_id[index] = mol_start_idx + i + 1;

        // Randomizing positions of the rest beads within the chain spatially
        for (int j = 1; j < natom_perchain; ++j)
        {
            // from index 1(=2nd bead of each chain)
            index = start_idx + i * natom_perchain + j;
            // assign the atom type
            particle_type[index] = itype;
            // assign the molecule ID
            molecule_id[index] = mol_start_idx + i + 1;
            // create an entry for a bond
            bonds[ibond][0] = index - 1;
            bonds[ibond][1] = index;
            ibond += 1;
            // Add the bond info to both atom[index] and atom[index-1]
            atom_bond_list[index][0] = index - 1;
            atom_bond_list[index - 1][1] = index;
            // Randomly get theta angle(xy plane)
            theta = ran1(idum) * 2.0 * M_PI;
            // randomly get cos() value ranging from -1 to 1 to set 2Pi rotation angle from z axis
            phi = acos(2.0 * ran1(idum) - 1.0);
            // move this much in x direction relative to previous bead
            coords[index][0] = coords[index - 1][0] + b * cos(theta) * sin(phi);
            // move this much in y direction relative to previous bead
            coords[index][1] = coords[index - 1][1] + b * sin(theta) * sin(phi);
            // move this much in z direction relative to previous bead
            coords[index][2] = coords[index - 1][2] + b * cos(phi);

            // If our new chain goes outside of the periodic box, we wrap around to the other side of the box
            coords[index][0] = periodic(coords[index][0], box_size);
            coords[index][1] = periodic(coords[index][1], box_size);
            coords[index][2] = periodic(coords[index][2], box_size);
            // For debug:
            // printf("%d\t", index);
            // printf("%d\t", particle_type[index]);
            // printf("%d\t", molecule_id[index]);
            // printf("%lf\t", coords[index][0]);
            // printf("%lf\t", coords[index][1]);
            // printf("%lf\t", coords[index][2]);
            // printf("\n");
            // printf("%d %d\n", bonds[ibond][0], bonds[ibond][1]);

        } // iterate to complete one chain
    }     // done initialization
}

// initialization of bottlebrush polymers
void initialize_bottlebrush(long long int nchain, int itype, double bb, int Nbb, int Nsc, int f_branch, double box_size,
                            int start_idx, int mol_start_idx, int *particle_type, int *molecule_id,
                            double **coords, int **bonds, int **atom_bond_list, long *idum)
{
    double theta, phi;
    int index, start_thischain;
    double b = sqrt(bb);
    int ibond;
    int natom_perchain;
    ibond = start_idx - mol_start_idx;
    natom_perchain = Nbb * (f_branch * Nsc + 1);
    for (int i = 0; i < nchain; ++i)
    {
        // Randomizing positions/locations of the first bead(=index 0, 1N, 2N-th...) of each polymer
        start_thischain = start_idx + i * natom_perchain;
        index = start_thischain;
        coords[index][0] = ran1(idum) * box_size;
        coords[index][1] = ran1(idum) * box_size;
        coords[index][2] = ran1(idum) * box_size;
        particle_type[index] = itype;
        molecule_id[index] = mol_start_idx + i + 1;

        // Create the backbone first
        for (int j = 1; j < Nbb; ++j)
        {
            // from index 1(=2nd bead of each chain)
            index = start_thischain + j;
            // assign the atom type
            particle_type[index] = itype;
            // assign the molecule ID
            molecule_id[index] = mol_start_idx + i + 1;
            // create an entry for a bond
            bonds[ibond][0] = index - 1;
            bonds[ibond][1] = index;
            ibond += 1;
            // Add the bond info to both atom[index] and atom[index-1]
            atom_bond_list[index][0] = index - 1;
            atom_bond_list[index - 1][1] = index;
            // Randomly get theta angle(xy plane)
            theta = ran1(idum) * 2.0 * M_PI;
            // randomly get cos() value ranging from -1 to 1 to set 2Pi rotation angle from z axis
            phi = acos(2.0 * ran1(idum) - 1.0);
            // move this much in x direction relative to previous bead
            coords[index][0] = coords[index - 1][0] + b * cos(theta) * sin(phi);
            // move this much in y direction relative to previous bead
            coords[index][1] = coords[index - 1][1] + b * sin(theta) * sin(phi);
            // move this much in z direction relative to previous bead
            coords[index][2] = coords[index - 1][2] + b * cos(phi);

            // If our new chain goes outside of the periodic box, we wrap around to the other side of the box
            coords[index][0] = periodic(coords[index][0], box_size);
            coords[index][1] = periodic(coords[index][1], box_size);
            coords[index][2] = periodic(coords[index][2], box_size);
        }

        // then for each backbone bead, attach f side-chains to it
        for (int ibb = 0; ibb < Nbb; ++ibb)
        {
            for (int i_branch = 0; i_branch < f_branch; ++i_branch)
            {
                index += 1;
                particle_type[index] = itype;
                molecule_id[index] = mol_start_idx + i + 1;
                bonds[ibond][0] = start_thischain + ibb;
                bonds[ibond][1] = index;
                ibond += 1;
                atom_bond_list[index][0] = start_thischain + ibb;
                atom_bond_list[start_thischain + ibb][i_branch + 2] = index;
                theta = ran1(idum) * 2.0 * M_PI;
                phi = acos(2.0 * ran1(idum) - 1.0);
                coords[index][0] = coords[start_thischain + ibb][0] + b * cos(theta) * sin(phi);
                coords[index][1] = coords[start_thischain + ibb][1] + b * sin(theta) * sin(phi);
                coords[index][2] = coords[start_thischain + ibb][2] + b * cos(phi);
                coords[index][0] = periodic(coords[index][0], box_size);
                coords[index][1] = periodic(coords[index][1], box_size);
                coords[index][2] = periodic(coords[index][2], box_size);

                for (int k = 1; k < Nsc; ++k)
                {
                    index += 1;
                    particle_type[index] = itype;
                    molecule_id[index] = mol_start_idx + i + 1;
                    bonds[ibond][0] = index - 1;
                    bonds[ibond][1] = index;
                    ibond += 1;
                    atom_bond_list[index][0] = index - 1;
                    atom_bond_list[index - 1][1] = index;
                    theta = ran1(idum) * 2.0 * M_PI;
                    phi = acos(2.0 * ran1(idum) - 1.0);
                    coords[index][0] = coords[index - 1][0] + b * cos(theta) * sin(phi);
                    coords[index][1] = coords[index - 1][1] + b * sin(theta) * sin(phi);
                    coords[index][2] = coords[index - 1][2] + b * cos(phi);
                    coords[index][0] = periodic(coords[index][0], box_size);
                    coords[index][1] = periodic(coords[index][1], box_size);
                    coords[index][2] = periodic(coords[index][2], box_size);
                }
            }
        }
    }
}
// // initialization of bottlebrush blockcopolymers
// void initialize_bbcp(long long int nchain, int itype, double bb, int Nbb_a, int Nbb_b, int Nsc_a, int Nsc_b, int f_brancha, int f_branchb, double box_size,
//                             int start_idx, int mol_start_idx, int *particle_typea, int *particle_typeb, int *molecule_id,
//                             double **coordsa, double **coordsb, int **bondsa, int **bondsb, int **atom_bond_lista, int **atom_bond_listb, long *idum)
// {
//     double theta, phi;
//     int index, start_thischain;
//     double b = sqrt(bb);
//     int ibond;
//     int natom_perchaina;
//     int natom_perchainb;
//     ibond = start_idx - mol_start_idx;
//     natom_perchaina = Nbb_a * (f_brancha * Nsc_a + 1);
//     natom_perchainb = Nbb_b * (f_branchb * Nsc_b + 1);
//     for (int i = 0; i < nchain; ++i)
//     {
//         // Randomizing positions/locations of the first bead(=index 0, 1N, 2N-th...) of each polymer
//         start_thischain = start_idx + i * natom_perchaina;
//         index = start_thischain;
//         coords[index][0] = ran1(idum) * box_size;
//         coords[index][1] = ran1(idum) * box_size;
//         coords[index][2] = ran1(idum) * box_size;
//         particle_type[index] = itype;
//         molecule_id[index] = mol_start_idx + i + 1;

//         // Create the backbone first
//         for (int j = 1; j < Nbb_a; ++j)
//         {
//             // from index 1(=2nd bead of each chain)
//             index = start_thischain + j;
//             // assign the atom type
//             particle_type[index] = itype;
//             // assign the molecule ID
//             molecule_id[index] = mol_start_idx + i + 1;
//             // create an entry for a bond
//             bonds[ibond][0] = index - 1;
//             bonds[ibond][1] = index;
//             ibond += 1;
//             // Add the bond info to both atom[index] and atom[index-1]
//             atom_bond_list[index][0] = index - 1;
//             atom_bond_list[index - 1][1] = index;
//             // Randomly get theta angle(xy plane)
//             theta = ran1(idum) * 2.0 * M_PI;
//             // randomly get cos() value ranging from -1 to 1 to set 2Pi rotation angle from z axis
//             phi = acos(2.0 * ran1(idum) - 1.0);
//             // move this much in x direction relative to previous bead
//             coords[index][0] = coords[index - 1][0] + b * cos(theta) * sin(phi);
//             // move this much in y direction relative to previous bead
//             coords[index][1] = coords[index - 1][1] + b * sin(theta) * sin(phi);
//             // move this much in z direction relative to previous bead
//             coords[index][2] = coords[index - 1][2] + b * cos(phi);

//             // If our new chain goes outside of the periodic box, we wrap around to the other side of the box
//             coords[index][0] = periodic(coords[index][0], box_size);
//             coords[index][1] = periodic(coords[index][1], box_size);
//             coords[index][2] = periodic(coords[index][2], box_size);
//         }

//         // then for each backbone bead, attach f side-chains to it
//         for (int ibb = 0; ibb < Nbb_a; ++ibb)
//         {
//             for (int i_branch = 0; i_branch < f_brancha; ++i_branch)
//             {
//                 index += 1;
//                 particle_type[index] = itype;
//                 molecule_id[index] = mol_start_idx + i + 1;
//                 bonds[ibond][0] = start_thischain + ibb;
//                 bonds[ibond][1] = index;
//                 ibond += 1;
//                 atom_bond_list[index][0] = start_thischain + ibb;
//                 atom_bond_list[start_thischain + ibb][i_branch + 2] = index;
//                 theta = ran1(idum) * 2.0 * M_PI;
//                 phi = acos(2.0 * ran1(idum) - 1.0);
//                 coords[index][0] = coords[start_thischain + ibb][0] + b * cos(theta) * sin(phi);
//                 coords[index][1] = coords[start_thischain + ibb][1] + b * sin(theta) * sin(phi);
//                 coords[index][2] = coords[start_thischain + ibb][2] + b * cos(phi);
//                 coords[index][0] = periodic(coords[index][0], box_size);
//                 coords[index][1] = periodic(coords[index][1], box_size);
//                 coords[index][2] = periodic(coords[index][2], box_size);

//                 for (int k = 1; k < Nsc_a; ++k)
//                 {
//                     index += 1;
//                     particle_type[index] = itype;
//                     molecule_id[index] = mol_start_idx + i + 1;
//                     bonds[ibond][0] = index - 1;
//                     bonds[ibond][1] = index;
//                     ibond += 1;
//                     atom_bond_list[index][0] = index - 1;
//                     atom_bond_list[index - 1][1] = index;
//                     theta = ran1(idum) * 2.0 * M_PI;
//                     phi = acos(2.0 * ran1(idum) - 1.0);
//                     coords[index][0] = coords[index - 1][0] + b * cos(theta) * sin(phi);
//                     coords[index][1] = coords[index - 1][1] + b * sin(theta) * sin(phi);
//                     coords[index][2] = coords[index - 1][2] + b * cos(phi);
//                     coords[index][0] = periodic(coords[index][0], box_size);
//                     coords[index][1] = periodic(coords[index][1], box_size);
//                     coords[index][2] = periodic(coords[index][2], box_size);
//                 }
//             }
//         }

//         // Another block attatched 
//         for (int j = 1; j < Nbb_b; ++j)
//         {
//             // from index 1(=2nd bead of each chain)
//             index += 1;
//             // assign the atom type
//             particle_type[index] = itype+1;
//             // assign the molecule ID
//             molecule_id[index] = mol_start_idx + i + 1;
//             // create an entry for a bond
//             if(j==1)
//             {
//                 bonds[ibond][0] = start_thischain + Nbb_a;
//                 bonds[ibond][1] = index;
//             }
//             else
//             {
//                 bonds[ibond][0] = index - 1;
//                 bonds[ibond][1] = index;
//             }

//             ibond += 1;
//             // Add the bond info to both atom[index] and atom[index-1]
//             atom_bond_list[index][0] = index - 1;
//             atom_bond_list[index - 1][1] = index;
//             // Randomly get theta angle(xy plane)
//             theta = ran1(idum) * 2.0 * M_PI;
//             // randomly get cos() value ranging from -1 to 1 to set 2Pi rotation angle from z axis
//             phi = acos(2.0 * ran1(idum) - 1.0);
//             // move this much in x direction relative to previous bead
//             coords[index][0] = coords[index - 1][0] + b * cos(theta) * sin(phi);
//             // move this much in y direction relative to previous bead
//             coords[index][1] = coords[index - 1][1] + b * sin(theta) * sin(phi);
//             // move this much in z direction relative to previous bead
//             coords[index][2] = coords[index - 1][2] + b * cos(phi);

//             // If our new chain goes outside of the periodic box, we wrap around to the other side of the box
//             coords[index][0] = periodic(coords[index][0], box_size);
//             coords[index][1] = periodic(coords[index][1], box_size);
//             coords[index][2] = periodic(coords[index][2], box_size);
//         }
//         // then for each backbone bead, attach f side-chains to it
//         for (int ibb = 0; ibb < Nbb_b; ++ibb)
//         {
//             for (int i_branch = 0; i_branch < f_branchb; ++i_branch)
//             {
//                 index += 1;
//                 particle_type[index] = itype+1;
//                 molecule_id[index] = mol_start_idx + i + 1;
//                 bonds[ibond][0] = start_thischain + Nbb_a+ibb;
//                 bonds[ibond][1] = index;
//                 ibond += 1;
//                 atom_bond_list[index][0] = start_thischain + ibb;
//                 atom_bond_list[start_thischain + ibb][i_branch + 2] = index;
//                 theta = ran1(idum) * 2.0 * M_PI;
//                 phi = acos(2.0 * ran1(idum) - 1.0);
//                 coords[index][0] = coords[start_thischain + ibb][0] + b * cos(theta) * sin(phi);
//                 coords[index][1] = coords[start_thischain + ibb][1] + b * sin(theta) * sin(phi);
//                 coords[index][2] = coords[start_thischain + ibb][2] + b * cos(phi);
//                 coords[index][0] = periodic(coords[index][0], box_size);
//                 coords[index][1] = periodic(coords[index][1], box_size);
//                 coords[index][2] = periodic(coords[index][2], box_size);

//                 for (int k = 1; k < Nsc_b; ++k)
//                 {
//                     index += 1;
//                     particle_type[index] = itype+1;
//                     molecule_id[index] = mol_start_idx + i + 1;
//                     bonds[ibond][0] = index - 1;
//                     bonds[ibond][1] = index;
//                     ibond += 1;
//                     atom_bond_list[index][0] = index - 1;
//                     atom_bond_list[index - 1][1] = index;
//                     theta = ran1(idum) * 2.0 * M_PI;
//                     phi = acos(2.0 * ran1(idum) - 1.0);
//                     coords[index][0] = coords[index - 1][0] + b * cos(theta) * sin(phi);
//                     coords[index][1] = coords[index - 1][1] + b * sin(theta) * sin(phi);
//                     coords[index][2] = coords[index - 1][2] + b * cos(phi);
//                     coords[index][0] = periodic(coords[index][0], box_size);
//                     coords[index][1] = periodic(coords[index][1], box_size);
//                     coords[index][2] = periodic(coords[index][2], box_size);
//                 }
//             }
//         }

//     }
// }
void initialize_brushblend(long long int nchain1, long long int nchain2, int itype1, int itype2, double bb1, double bb2, int Nbb1, int Nbb2, int Nsc1, int Nsc2, int f_branch1, int f_branch2, double box_size,
                            int start_idx1, int mol_start_idx1, int *particle_type, int *bond_type, int *molecule_id,
                            double **coords, int **bonds, int **atom_bond_list, long *idum)
{
    double theta, phi;
    int index, start_thischain;
    int start_idx2, mol_start_idx2;
    double b1 = sqrt(bb1);
    double b2 = sqrt(bb2);
    int ibond;
    int natom_perchain1, natom_perchain2;
    ibond = start_idx1 - mol_start_idx1;
    natom_perchain1 = Nbb1 * (f_branch1 * Nsc1 + 1);
    for (int i = 0; i < nchain1; ++i)
    {
        // Randomizing positions/locations of the first bead(=index 0, 1N, 2N-th...) of each polymer
        start_thischain = start_idx1 + i * natom_perchain1;
        index = start_thischain;
        coords[index][0] = ran1(idum) * box_size;
        coords[index][1] = ran1(idum) * box_size;
        coords[index][2] = ran1(idum) * box_size;
        particle_type[index] = itype1;
        molecule_id[index] = mol_start_idx1 + i + 1;

        // Create the type 1 backbone first
        for (int j = 1; j < Nbb1; ++j)
        {
            // from index 1(=2nd bead of each chain)
            index = start_thischain + j;
            particle_type[index] = itype1;
            bond_type[ibond] = itype1+1;
            molecule_id[index] = mol_start_idx1 + i + 1;
            bonds[ibond][0] = index - 1;
            bonds[ibond][1] = index;
            ibond += 1;
            atom_bond_list[index][0] = index - 1;
            atom_bond_list[index - 1][1] = index;
            theta = ran1(idum) * 2.0 * M_PI;
            phi = acos(2.0 * ran1(idum) - 1.0);
            coords[index][0] = coords[index - 1][0] + b1 * cos(theta) * sin(phi);
            coords[index][1] = coords[index - 1][1] + b1 * sin(theta) * sin(phi);
            coords[index][2] = coords[index - 1][2] + b1 * cos(phi);

            coords[index][0] = periodic(coords[index][0], box_size);
            coords[index][1] = periodic(coords[index][1], box_size);
            coords[index][2] = periodic(coords[index][2], box_size);
        }

        // then for each backbone bead, attach f side-chains to it
        for (int ibb = 0; ibb < Nbb1; ++ibb)
        {
            for (int i_branch = 0; i_branch < f_branch1; ++i_branch)
            {
                index += 1;
                particle_type[index] = itype1;
                bond_type[ibond] = itype1+1;
                molecule_id[index] = mol_start_idx1 + i + 1;
                bonds[ibond][0] = start_thischain + ibb;
                bonds[ibond][1] = index;
                ibond += 1;
                atom_bond_list[index][0] = start_thischain + ibb;
                atom_bond_list[start_thischain + ibb][i_branch + 2] = index;
                theta = ran1(idum) * 2.0 * M_PI;
                phi = acos(2.0 * ran1(idum) - 1.0);
                coords[index][0] = coords[start_thischain + ibb][0] + b1 * cos(theta) * sin(phi);
                coords[index][1] = coords[start_thischain + ibb][1] + b1 * sin(theta) * sin(phi);
                coords[index][2] = coords[start_thischain + ibb][2] + b1 * cos(phi);
                coords[index][0] = periodic(coords[index][0], box_size);
                coords[index][1] = periodic(coords[index][1], box_size);
                coords[index][2] = periodic(coords[index][2], box_size);

                for (int k = 1; k < Nsc1; ++k)
                {
                    index += 1;
                    particle_type[index] = itype1;
                    bond_type[ibond] = itype1+1;
                    molecule_id[index] = mol_start_idx1 + i + 1;
                    bonds[ibond][0] = index - 1;
                    bonds[ibond][1] = index;
                    ibond += 1;
                    atom_bond_list[index][0] = index - 1;
                    atom_bond_list[index - 1][1] = index;
                    theta = ran1(idum) * 2.0 * M_PI;
                    phi = acos(2.0 * ran1(idum) - 1.0);
                    coords[index][0] = coords[index - 1][0] + b1 * cos(theta) * sin(phi);
                    coords[index][1] = coords[index - 1][1] + b1 * sin(theta) * sin(phi);
                    coords[index][2] = coords[index - 1][2] + b1 * cos(phi);
                    coords[index][0] = periodic(coords[index][0], box_size);
                    coords[index][1] = periodic(coords[index][1], box_size);
                    coords[index][2] = periodic(coords[index][2], box_size);
                }
            }
        }
    }

    // Now initialize the type 2 backbone
    natom_perchain2 = Nbb2 * (f_branch2 * Nsc2 + 1);
    start_idx2 = start_idx1 + nchain1 * natom_perchain1;
    mol_start_idx2 = mol_start_idx1 + nchain1;
   for (int i = 0; i < nchain2; ++i)
    {
        // Randomizing positions/locations of the first bead(=index 0, 1N, 2N-th...) of each polymer
        start_thischain = start_idx2 + i * natom_perchain2;
        index = start_thischain;
        coords[index][0] = ran1(idum) * box_size;
        coords[index][1] = ran1(idum) * box_size;
        coords[index][2] = ran1(idum) * box_size;
        particle_type[index] = itype2;
        molecule_id[index] = mol_start_idx2 + i + 1;

        // Create the type 1 backbone first
        for (int j = 1; j < Nbb2; ++j)
        {
            index = start_thischain + j;
            particle_type[index] = itype2;
            bond_type[ibond] = itype2+1;
            molecule_id[index] = mol_start_idx2 + i + 1;
            bonds[ibond][0] = index - 1;
            bonds[ibond][1] = index;
            ibond += 1;
            atom_bond_list[index][0] = index - 1;
            atom_bond_list[index - 1][1] = index;
            theta = ran1(idum) * 2.0 * M_PI;
            phi = acos(2.0 * ran1(idum) - 1.0);
            coords[index][0] = coords[index - 1][0] + b2 * cos(theta) * sin(phi);
            coords[index][1] = coords[index - 1][1] + b2 * sin(theta) * sin(phi);
            coords[index][2] = coords[index - 1][2] + b2 * cos(phi);

            coords[index][0] = periodic(coords[index][0], box_size);
            coords[index][1] = periodic(coords[index][1], box_size);
            coords[index][2] = periodic(coords[index][2], box_size);
        }

        // then for each backbone bead, attach f side-chains to it
        for (int ibb = 0; ibb < Nbb2; ++ibb)
        {
            for (int i_branch = 0; i_branch < f_branch2; ++i_branch)
            {
                index += 1;
                particle_type[index] = itype2;
                bond_type[ibond] = itype2+1;
                molecule_id[index] = mol_start_idx2 + i + 1;
                bonds[ibond][0] = start_thischain + ibb;
                bonds[ibond][1] = index;
                ibond += 1;
                atom_bond_list[index][0] = start_thischain + ibb;
                atom_bond_list[start_thischain + ibb][i_branch + 2] = index;
                theta = ran1(idum) * 2.0 * M_PI;
                phi = acos(2.0 * ran1(idum) - 1.0);
                coords[index][0] = coords[start_thischain + ibb][0] + b2 * cos(theta) * sin(phi);
                coords[index][1] = coords[start_thischain + ibb][1] + b2 * sin(theta) * sin(phi);
                coords[index][2] = coords[start_thischain + ibb][2] + b2 * cos(phi);
                coords[index][0] = periodic(coords[index][0], box_size);
                coords[index][1] = periodic(coords[index][1], box_size);
                coords[index][2] = periodic(coords[index][2], box_size);

                for (int k = 1; k < Nsc2; ++k)
                {
                    index += 1;
                    particle_type[index] = itype2;
                    bond_type[ibond] = itype2+1;
                    molecule_id[index] = mol_start_idx2 + i + 1;
                    bonds[ibond][0] = index - 1;
                    bonds[ibond][1] = index;
                    ibond += 1;
                    atom_bond_list[index][0] = index - 1;
                    atom_bond_list[index - 1][1] = index;
                    theta = ran1(idum) * 2.0 * M_PI;
                    phi = acos(2.0 * ran1(idum) - 1.0);
                    coords[index][0] = coords[index - 1][0] + b2 * cos(theta) * sin(phi);
                    coords[index][1] = coords[index - 1][1] + b2 * sin(theta) * sin(phi);
                    coords[index][2] = coords[index - 1][2] + b2 * cos(phi);
                    coords[index][0] = periodic(coords[index][0], box_size);
                    coords[index][1] = periodic(coords[index][1], box_size);
                    coords[index][2] = periodic(coords[index][2], box_size);
                }
            }
        }
    }
            
}

// int main(int argc, const char *argv[])
// {
//     long long int natom_perchain, natom_total, natom_A, natom_B, Nbb, Nsc, f_branch;
//     long long int nchain, nchain_A, nchain_B;
//     int nbonds_perchain, nbonds_total, nint;
//     int Nbb1, Nbb2, Nsc1, Nsc2, f_branch1, f_branch2;
//     int nchain1, nchain2;

//     int pm_type;
//     double kappa, kspring, chi, bondlen, bondlen_sqr, Nbar;
//     double f_A;
//     double grid_size;
//     int maxsite_1d;

//     int save_every, tot_mc_cycles, equil_cycles, grid_shift_every;
//     double ave_bondlen, reedsq_bb, reedsq_sc, rgsq_bb, rgsq_sc;
//     time_t start_t, current_t;
//     double diff_t;
//     int diff_days, diff_hours, diff_minutes, diff_seconds;
//     time(&start_t);
//     Nbb1=20;
//     Nbb2=20;
//     Nsc1=3;
//     Nsc2=3;
//     f_branch1=1;
//     f_branch2=1;
//     nchain1=1;
//     nchain2=1;
//     natom_total = nchain1*(Nbb1*(f_branch1*Nsc1+1))+nchain2*(Nbb2*(f_branch2*Nsc2+1));
//     nbonds_total = natom_total - nchain1 - nchain2;
//     // Memory allocation
//     // type of particle, 1d array (natoms)
//     int *particle_type = calloc(natom_total, sizeof(int));
//     int *bond_type = calloc(nbonds_total, sizeof(int));
//     // molecule ID, 1d array (natoms)
//     int *molecule_id = calloc(natom_total, sizeof(int));
//     // coordinates of particles, 2d array (natoms, 3)
//     double **coords;
//     coords = (double **)calloc(natom_total, sizeof(double *));
//     for (int i = 0; i < natom_total; ++i)
//     {
//         coords[i] = (double *)calloc(3, sizeof(double));
//     }
//     // coords[0][0] = 0;
//     // bond table, 2d array (nbonds, 3)
//     int **bonds;
//     bonds = (int **)calloc(nbonds_total, sizeof(int *));
//     for (int i = 0; i < nbonds_total; ++i)
//     {
//         bonds[i] = (int *)calloc(2, sizeof(int));
//     }
//     // atom bond list, each entry lists the atoms that the given atom is bonded with
//     int **atom_bond_list;
//     atom_bond_list = (int **)calloc(natom_total, sizeof(int *));
//     for (int i = 0; i < natom_total; ++i)
//     {
//         // Assuming maximum 10 bonds per atom!!!
//         atom_bond_list[i] = (int *)calloc(10, sizeof(int));
//     }
//     // Then initialize this list to -1 to make sure check can be performed
//     for (int i = 0; i < natom_total; ++i)
//     {
//         // Assuming maximum 10 bonds per atom!!!
//         for (int j = 0; j < 10; ++j)
//         {
//             atom_bond_list[i][j] = -1;
//         }
//     }
//     long *idum = malloc(sizeof(long));
//     *idum = initRan();
//     // Initialize polymers
//     int start_idx = 0;
//     int mol_start_idx = 0;
//     initialize_brushblend(nchain1, nchain2, 0, 1, 1, 1, Nbb1, Nbb2, Nsc1, Nsc2, f_branch1, f_branch2, 20,
//                           start_idx, mol_start_idx, particle_type, bond_type, molecule_id,
//                           coords, bonds, atom_bond_list, idum);

//     // Write the initial configuration to a file
//     FILE *init_config;
//     init_config = fopen("ini.cfg", "w");

//     fprintf(init_config, "#Initial configuration\n");
//     fprintf(init_config, "%lld atoms\n", natom_total);
//     fprintf(init_config, "2 atom types\n");
//     fprintf(init_config, "%d bonds\n", nbonds_total);
//     fprintf(init_config, "2 bond types\n");
//     fprintf(init_config, "0.0 %lf xlo xhi\n", 20.0);
//     fprintf(init_config, "0.0 %lf ylo yhi\n", 20.0);
//     fprintf(init_config, "0.0 %lf zlo zhi\n", 20.0);
//     fprintf(init_config, "0 0 0 xy xz yz\n");
//     fprintf(init_config, "\n");
//     fprintf(init_config, "Atoms # bond\n");
//     fprintf(init_config, "\n");

//     for (int i = 0; i < natom_total; ++i)
//     {
//         fprintf(init_config, "%d %d %d %f %f %f\n", i + 1, molecule_id[i], particle_type[i] + 1,
//                 coords[i][0], coords[i][1], coords[i][2]);
//     }

//     fprintf(init_config, "\n");
//     fprintf(init_config, "Bonds\n");
//     fprintf(init_config, "\n");

//     for (int i = 0; i < nbonds_total; ++i)
//     {
//         fprintf(init_config, "%d %d %d %d\n", i + 1, bond_type[i], bonds[i][0] + 1, bonds[i][1] + 1);
//     }
//     fclose(init_config);
// }

// float ran1(long *idum)
// {
//     int j;
//     long k;
//     static long idum2 = 123456789;
//     static long iy = 0;
//     static long iv[NTAB];
//     float temp;

//     if (*idum <= 0)
//     {
//         if (-(*idum) < 1)
//             *idum = 1;
//         else
//             *idum = -(*idum);
//         idum2 = (*idum);
//         for (j = NTAB + 7; j >= 0; --j)
//         {
//             k = (*idum) / IQ1;
//             *idum = IA1 * (*idum - k * IQ1) - k * IR1;
//             if (*idum < 0)
//                 *idum += IM1;
//             if (j < NTAB)
//                 iv[j] = *idum;
//         }
//         iy = iv[0];
//     }
//     k = (*idum) / IQ1;
//     *idum = IA1 * (*idum - k * IQ1) - k * IR1;
//     if (*idum < 0)
//         *idum += IM1;
//     k = idum2 / IQ2;
//     if (*idum < 0)
//         idum2 += IM2;
//     j = iy / NDIV;
//     iy = iv[j] - idum2;
//     iv[j] = *idum;
//     if (iy < 1)
//         iy += IMM1;
//     if ((temp = AM * iy) > RNMX)
//         return RNMX;
//     else
//         return temp;
// }

// // random number with normal distribution N(0,1)
// float gasdev(long *idum)
// {
//     float ran1(long *idum);
//     static int iset = 0;
//     static float gset;
//     float fac, rsq, v1, v2;

//     if (*idum < 0)
//         iset = 0;
//     if (iset == 0)
//     {
//         do
//         {
//             v1 = 2.0 * ran1(idum) - 1.0;
//             v2 = 2.0 * ran1(idum) - 1.0;
//             rsq = v1 * v1 + v2 * v2;
//         } while (rsq >= 1.0 || rsq == 0.0);
//         fac = sqrt(-2.0 * log(rsq) / rsq);
//         gset = v1 * fac;
//         iset = 1;
//         return v2 * fac;
//     }
//     else
//     {
//         iset = 0;
//         return gset;
//     }
// }

// // initialize random number seed
// long initRan()
// {
//     time_t seconds;
//     time(&seconds);
//     return -1 * (unsigned long)(seconds);
// }

// double periodic(double ri, double bl)
// {
//     if (ri < 0)
//         ri += bl;
//     else if (ri > bl)
//         ri -= bl;

//     return ri;
// }