import pandas as pd
import numpy as np
import sys
def parse_lammps_dump_all_timesteps(filename):
    # Store positions for each timestep
    timesteps = []
    current_positions = []
    with open(filename) as f:
        lines = f.readlines()
    
    in_snapshot = False

    for line in lines:
        if "ITEM: TIMESTEP" in line:
            if current_positions:
                # Store the positions of the last timestep before starting a new one
                timesteps.append(pd.DataFrame(current_positions, columns=['id', 'mol', 'type', 'x', 'y', 'z']))
                current_positions = []
                in_snapshot = False
        elif "ITEM: ATOMS" in line:
            in_snapshot = True
        elif "ITEM:" in line:
            in_snapshot = False
        elif in_snapshot:
            # Adjust based on your dump file format (assuming columns are id, x, y, z, etc.)
            parts = line.strip().split()
            # print(parts)
            current_positions.append([int(parts[0]), int(parts[1]), int(parts[2]), float(parts[3]), float(parts[4]), float(parts[5])])
    
    # Append the last timestep
    if current_positions:
        timesteps.append(pd.DataFrame(current_positions, columns=['id', 'mol', 'type', 'x', 'y', 'z']))

    return timesteps
## Test parse_lammps_dump_all_timesteps
timesteps = parse_lammps_dump_all_timesteps(sys.argv[1])
print(timesteps[-1])

def structure_factor(positions, q_values, lx, ly, lz):
    num_atoms = len(positions)
    s_q = []

    # Calculate the structure factor for each q in q_values
    for q in q_values:
        real_sum = 0.0
        imag_sum = 0.0
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                # Distance vector between atoms i and j
                r_ij = positions.iloc[i] - positions.iloc[j]
                # print(r_ij)
                # Apply minimum image convention for each component
                r_ij[0] -= np.rint(r_ij[0] / lx) * lx  # x component
                r_ij[1] -= np.rint(r_ij[1] / ly) * ly  # y component
                r_ij[2] -= np.rint(r_ij[2] / lz) * lz  # z component
                
                # Calculate distance magnitude
                r = np.linalg.norm(r_ij)
                # print(r)
                # Calculate real and imaginary components of the Fourier transform
                real_sum += np.cos(q * r)
                imag_sum += np.sin(q * r)
        print(real_sum)
        print(imag_sum)
        # Normalize by the number of atoms and calculate S(q)
        s_q_value = 1.0 + (2.0 / num_atoms) * (real_sum**2 + imag_sum**2)**0.5
        s_q.append(s_q_value)
        print(s_q_value)
    
    return np.array(s_q)


q_values = np.arange(0.01, 100 , 20)

positions = timesteps[-1][['x', 'y', 'z']]
s_q = structure_factor(positions, q_values, 100, 100, 20)


