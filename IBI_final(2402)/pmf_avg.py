import sys
import numpy as np

# Specify the path to your text file and the desired iteration number
file_path = sys.argv[1]
desired_iteration = int(sys.argv[2])  # Take the iteration number from command line

# Read the file and split it into iterations
with open(file_path, 'r') as file:
    data = file.read().split('Iteration')

# Ensure the desired iteration is within the range of available iterations
if desired_iteration < 1 or desired_iteration >= len(data):
    print(f"Error: Desired iteration {desired_iteration} is out of range. Available iterations are from 1 to {len(data)-1}.")
    sys.exit(1)

# Get the data from the specified iteration
iteration_data = data[desired_iteration].strip().split('\n')[1:]

# Convert the string data into a NumPy array
iteration_array = np.array([list(map(float, line.split())) for line in iteration_data])

# Save the NumPy array to a text file
output_filename = f"finalpmf_iteration_{desired_iteration}.txt"
np.savetxt(output_filename, iteration_array, delimiter='\t', fmt="%.4f")