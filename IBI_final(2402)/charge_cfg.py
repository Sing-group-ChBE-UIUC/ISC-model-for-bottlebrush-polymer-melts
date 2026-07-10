import sys
# Input and output file names
input_file = sys.argv[1]  # Replace with your input file name
output_file = 'charged.cfg'  # Replace with your desired output file name

# Read the input file
with open(input_file, 'r') as file:
    lines = file.readlines()

# Locate the start of the "Atoms" section
atom_start_index = lines.index('Atoms  # full\n') + 2  # Start after the 'Atoms' header

# Find where the "Atoms" section ends and the "Bonds" section begins
bonds_start_index = None
for i in range(atom_start_index, len(lines)):
    if lines[i].startswith('Bonds'):
        bonds_start_index = i
        break

# Get only the atoms section lines
atoms_section = lines[atom_start_index:bonds_start_index]

# Parse the atoms information
atoms_info = []
for line in atoms_section:
    if line.strip():  # Skip empty lines
        parts = line.split()
        atom_id = int(parts[0])
        molecule_id = int(parts[1])
        atom_type = int(parts[2])
        charge = float(parts[3])
        x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
        atoms_info.append([atom_id, molecule_id, atom_type, charge, x, y, z])

# Sort atoms by molecule ID and then by atom ID
atoms_info.sort(key=lambda x: (x[1], x[0]))

# Assign new charges
last_molecule_id = -1
first_atom_index = None

for i, atom in enumerate(atoms_info):
    molecule_id = atom[1]

    if molecule_id != last_molecule_id:
        # Assign charge 0.1 to the first monomer of the new molecule
        if first_atom_index is not None:
            atoms_info[first_atom_index][3] = 0.1  # Set charge of the first atom of the previous molecule
        first_atom_index = i

        last_molecule_id = molecule_id

    # Set charge of the last monomer of the molecule to -0.1
    if i == len(atoms_info) - 1 or atoms_info[i + 1][1] != molecule_id:
        atom[3] = -0.1  # Last monomer of the molecule
    else:
        atom[3] = 0.0  # Reset other charges

# Write the updated data to a new file
with open(output_file, 'w') as file:
    # Write the header part
    file.writelines(lines[:atom_start_index])

    # Write the updated atoms section
    for atom in atoms_info:
        file.write(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]:.1f} {atom[4]} {atom[5]} {atom[6]}\n")

    # Write the rest of the file (starting from the Bonds section)
    file.write('\n')
    file.writelines(lines[bonds_start_index:])