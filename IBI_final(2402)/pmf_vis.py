import numpy as np
import matplotlib.pyplot as plt
import sys

# Read the data from a file (replace 'data.txt' with your file name)
with open(sys.argv[1], 'r') as file:
    lines = file.readlines()

# Initialize lists to store the data and a variable to track the iteration number
data = []
iteration_number = 0
frame = int(sys.argv[3])

# Iterate through the lines of the file to extract the data
for line in lines:
    if line.startswith('Iteration'):
        # A new iteration starts
        if data:
            if iteration_number <frame and iteration_number >0 :
            #if iteration_number == int(sys.argv[3]):
                distances, gr_values, mf = zip(*data)
                #
                plt.plot(distances, gr_values)
        data = []
        iteration_number += 1
    else:
        # Split the line into distance and g(r) values
        parts = line.split()
        if len(parts) == 3:
            distance, gr, mf = map(float, parts)
            data.append((distance, gr, mf))

# Append the last iteration's data
if data:
    if iteration_number <frame and iteration_number >0:
        distances, gr_values, mf = zip(*data)
        plt.plot(distances, gr_values)

index, r, gro, mf = np.loadtxt(sys.argv[2], unpack=True)
plt.plot(r, gro, label="target",linewidth=2.5, color='black')
plt.xlabel('Distance')
plt.ylabel('V(r)')
plt.ylim(-0.5, 8)
#plt.ylabel('mf')
#plt.title('Plot for All Iterations (excluding first iteration)')
#plt.legend()
#plt.xlim(0,13)
plt.savefig("pmf.pdf")
plt.show()
plt.close()
