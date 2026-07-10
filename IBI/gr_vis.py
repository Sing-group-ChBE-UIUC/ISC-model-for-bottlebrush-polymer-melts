import numpy as np
import matplotlib.pyplot as plt
import sys

# font = {'family': 'serif', 'serif': ['Computer Modern'], 'size': 13}
# plt.rc('font', **font)
# plt.rc('text', usetex=True)
# plt.rcParams['text.usetex'] = False
cmap1 = plt.get_cmap('Blues')
cmap2 = plt.get_cmap('Oranges')
# Read the data from a file (replace 'data.txt' with your file name)
# with open(sys.argv[1], 'r') as file:
#     lines = file.readlines()

# # Initialize lists to store the data and a variable to track the iteration number
data = []
iteration_number = 0
frame = int(sys.argv[5])

# # Iterate through the lines of the file to extract the data
# for line in lines:
#     if line.startswith('Iteration'):
#         # A new iteration starts
#         if data:
#             print(iteration_number)
#             if iteration_number <frame and iteration_number > int(sys.argv[6]):
#             #if iteration_number == int(sys.argv[3]):
#                 distances, gr_values = zip(*data)
#                 #
#                 plt.plot(distances, gr_values, alpha=0.5, label="i=%d"%iteration_number)#color='blue'
#                 print(gr_values)
#         data = []
#         iteration_number += 1
#     else:
#         # Split the line into distance and g(r) values
#         parts = line.split()
#         if len(parts) == 3:
#             distance, gr = map(float, parts)
#             data.append((distance, gr))
#             print(gr)
# # Append the last iteration's data
# if data:
#     if iteration_number <frame and iteration_number >=int(sys.argv[6]):
#         distances, gr_values = zip(*data)
#         # plt.plot(distances, gr_values, color='blue', alpha=0.5, label=r"$g_{i,AA}(r)$")
#         # plt.plot(distances, grBB_values, label="BB")

#         r, gro= np.loadtxt(sys.argv[2], unpack=True)
#         plt.plot(r, gro, label=r"$g_{AA}(r)$",linewidth=2.5, color='blue')

with open(sys.argv[3], 'r') as file:
    lines = file.readlines()

# Initialize lists to store the data and a variable to track the iteration number
data = []
iteration_number = 0

# Iterate through the lines of the file to extract the data
for line in lines:
    if line.startswith('Iteration'):
        # A new iteration starts
        if data:
            if iteration_number <frame and iteration_number >=int(sys.argv[6]) :
            #if iteration_number == int(sys.argv[3]):
                distances, gr_values = zip(*data)
                # print(distances)
                plt.plot(distances, gr_values, alpha=0.5, label="i=%d"%iteration_number)#color='orange'

                # plt.plot(distances, gr_values, color=cmap2(0.1+float(iteration_number)/4))
        data = []
        iteration_number += 1
    else:
        # Split the line into distance and g(r) values
        parts = line.split()
        if len(parts) == 2:
            distance, gr = map(float, parts)
            data.append((distance, gr))

# # Append the last iteration's data
# if data:
#     if iteration_number <frame and iteration_number >=int(sys.argv[6]):
#         distances, gr_values = zip(*data)
#         plt.plot(distances, gr_values, color='orange', alpha=0.5, label=r"$g_{i,AB}(r)$")

r, gro= np.loadtxt(sys.argv[4], unpack=True)
plt.plot(r, gro, label=r"$g_{AB}(r)$",linewidth=2.5, color='orange')


plt.xlabel('r[b]')
plt.ylabel('g(r)')
plt.ylim(0, 3)
plt.xlim(0, 23)
#plt.ylabel('mf')
#plt.title('Plot for All Iterations (excluding first iteration)')
plt.legend(loc='lower right')
#plt.xlim(0,13)
plt.savefig("gr.pdf")
plt.show()
plt.close()
