import numpy as np
import matplotlib.pyplot as plt

def parse_lammps_dump(filename):
    """
    Parse LAMMPS dump file with multiple time steps, triclinic boxes, and atomic data.
    Returns:
        - positions_list: List of atomic positions for all time steps
        - types_list: List of particle types for all time steps
        - box_bounds_list: List of box bounds for all time steps
    """
    positions_list = []
    types_list = []
    box_bounds_list = []

    with open(filename, "r") as file:
        lines = file.readlines()

    i = 0
    while i < len(lines):
        if "ITEM: TIMESTEP" in lines[i]:
            i += 1  # Skip timestep value
            i += 1  # Skip "ITEM: NUMBER OF ATOMS"
            i += 1
            n_atoms = int(lines[i].strip())
            i += 1  # Skip "ITEM: BOX BOUNDS xy xz yz pp pp pp"
            i += 1
            # Parse box bounds (triclinic box handling)
            box_bounds = np.zeros((3, 2))
            tilt_factors = np.zeros(3)
            for j in range(3):
                box_line = list(map(float, lines[i].split()))
                box_bounds[j, :] = box_line[:2]
                if len(box_line) == 3:
                    tilt_factors[j] = box_line[2]
                i += 1

            # Parse atomic positions and types
            i += 1  # Skip "ITEM: ATOMS id mol type x y z"
            positions = []
            types = []
            for j in range(n_atoms):
                atom_data = list(map(float, lines[i].split()))
                positions.append(atom_data[3:6])  # Extract x, y, z
                types.append(int(atom_data[2]))   # Extract type
                i += 1

            positions_list.append(np.array(positions))
            types_list.append(np.array(types))
            box_bounds_list.append((box_bounds, tilt_factors))
        else:
            i += 1

    return positions_list, types_list, box_bounds_list

def extract_timestep(filename, target_timestep, output_file):
    """
    Extract data for a specific timestep from a LAMMPS dump file.

    Parameters:
        filename (str): Path to the input dump file.
        target_timestep (int): The timestep to extract.
        output_file (str): Path to the output file to save the extracted timestep.
    """
    with open(filename, "r") as file:
        lines = file.readlines()

    extracted_lines = []
    inside_target_timestep = False

    i = 0
    while i < len(lines):
        if "ITEM: TIMESTEP" in lines[i]:
            timestep = int(lines[i + 1].strip())
            if timestep == target_timestep:
                inside_target_timestep = True
                extracted_lines.append(lines[i])  # ITEM: TIMESTEP
                extracted_lines.append(lines[i + 1])  # Timestep value
            else:
                inside_target_timestep = False

        if inside_target_timestep:
            extracted_lines.append(lines[i])

        # Break out once the timestep data is fully extracted
        if inside_target_timestep and "ITEM: ATOMS" in lines[i]:
            for j in range(int(lines[i + 1])):
                extracted_lines.append(lines[i + 1 + j])
            break

        i += 1

    if extracted_lines:
        with open(output_file, "w") as out_file:
            out_file.writelines(extracted_lines)
        print(f"Timestep {target_timestep} extracted and saved to {output_file}")
    else:
        print(f"Timestep {target_timestep} not found in {filename}")
        
def convert_to_orthogonal(positions, box_bounds, tilt_factors):
    """
    Convert atomic positions from triclinic to orthogonal coordinates.
    """
    xlo, xhi = box_bounds[0]
    ylo, yhi = box_bounds[1]
    zlo, zhi = box_bounds[2]
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo

    # Construct triclinic transformation matrix
    h_matrix = np.array([
        [lx, tilt_factors[0], tilt_factors[1]],
        [0,  ly, tilt_factors[2]],
        [0,   0,  lz]
    ])

    # Normalize positions by the box dimensions
    fractional_positions = (positions - np.array([xlo, ylo, zlo])) / np.array([lx, ly, lz])
    orthogonal_positions = fractional_positions @ h_matrix
    return orthogonal_positions


def find_lamellar_direction(positions):
    """
    Find the lamellar direction from the atomic positions using principal component analysis (PCA).
    """
    cov_matrix = np.cov(positions, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
    lamellar_direction = eigenvectors[:, np.argmin(eigenvalues)]
    return lamellar_direction


def project_positions(positions, direction):
    """
    Project atomic positions onto the given direction vector.
    """
    direction = direction / np.linalg.norm(direction)
    projection = positions @ direction
    return projection


def calculate_box_length(box_bounds, tilt_factors, lamellar_direction):
    """
    Calculate the box length along the lamellar direction, considering tilt factors.
    """
    xlo, xhi = box_bounds[0]
    ylo, yhi = box_bounds[1]
    zlo, zhi = box_bounds[2]
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo

    # Construct the triclinic box vectors
    triclinic_vectors = np.array([
        [lx, 0, 0],
        [tilt_factors[0], ly, 0],
        [tilt_factors[1], tilt_factors[2], lz]
    ])

    # Calculate the projection of the lamellar direction onto the triclinic box vectors
    lamellar_direction_normalized = lamellar_direction / np.linalg.norm(lamellar_direction)
    box_length = np.linalg.norm(triclinic_vectors @ lamellar_direction_normalized)

    return box_length


def calculate_density_profile(projection, num_bins, box_length):
    """
    Calculate the density profile along the projected lamellar direction.
    """
    bin_edges = np.linspace(0, box_length, num_bins + 1)
    hist, _ = np.histogram(projection, bins=bin_edges)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    density = hist / (box_length / num_bins)
    return bin_centers, density


# Main script
if __name__ == "__main__":
    # Input file and parameters
    # dump_file = "dump.lammpstrj"
    dump_file = "dump_final.cfg"
    num_bins = 100

    # Parse dump file
    positions_list, types_list, box_bounds_list = parse_lammps_dump(dump_file)

    # Analyze each frame
    density_profiles = []
    for positions, types, (box_bounds, tilt_factors) in zip(positions_list, types_list, box_bounds_list):
        # Filter positions for type 1 particles
        type_1_positions = positions[types == 1]

        # Step 1: Convert positions to orthogonal coordinates
        orthogonal_positions = convert_to_orthogonal(type_1_positions, box_bounds, tilt_factors)

        # Step 2: Find lamellar direction for this frame
        lamellar_direction = find_lamellar_direction(orthogonal_positions)
        print("Lamellar direction:", lamellar_direction)
        lamellar_direction = np.array([1, 1, 0])  # Manually set lamellar direction

        # Step 3: Project positions onto the lamellar direction
        projection = project_positions(orthogonal_positions, lamellar_direction)

        # Step 4: Calculate the box length along the lamellar direction
        box_length = calculate_box_length(box_bounds, tilt_factors, lamellar_direction)
        box_length = 400

        # Step 5: Calculate the density profile for this frame
        bin_centers, density = calculate_density_profile(projection, num_bins, box_length)
        density_profiles.append(density)

    plt.plot(bin_centers, np.mean(density_profiles, axis=0))
    plt.show()
    
    # Save or visualize density profiles
    np.savetxt("density_profiles_type_1.txt", np.array(density_profiles))
    print("Density profiles for type 1 particles saved to density_profiles_type_1.txt")
