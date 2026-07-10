import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import skewnorm

cmap1 = plt.get_cmap('CMRmap')
# Initialize lists and arrays to store data
chiN = []
r1ab = []
chiN_data = {}
with open(sys.argv[1], "r") as file:
    # Step 1: Read the header line
    header = file.readline().strip().split()
    
    # The first column is 'r1ab', the rest are 'chiN' values
    chiN = [int(col.replace('chiN', '')) for col in header[1:]]
    
    header = file.readline().strip().split()
    header = file.readline().strip().split()
    header = file.readline().strip().split()
    # Step 2: Read the data lines
    data = np.loadtxt(file)
    
    # Separate the data into r1ab and chiN arrays
    r = data[:, 0]
    for i, chiN_value in enumerate(chiN):
        chiN_data[f'chiN{chiN_value}'] = data[:, i+1]

# Step 3: Assign chiN arrays to individual variables
hr1 = chiN_data[f'chiN{chiN[0]}']
hr2 = chiN_data[f'chiN{chiN[1]}']
hr3 = chiN_data[f'chiN{chiN[2]}']

# Initialize lists and arrays to store data
chiN = []
r1ab = []
chiN_data = {}
with open(sys.argv[2], "r") as file:
    # Step 1: Read the header line
    header = file.readline().strip().split()
    
    # The first column is 'r1ab', the rest are 'chiN' values
    chiN = [int(col.replace('chiN', '')) for col in header[1:]]
    
    # Step 2: Read the data lines
    data = np.loadtxt(file)
    
    # Separate the data into r1ab and chiN arrays
    r = data[:, 0]
    for i, chiN_value in enumerate(chiN):
        chiN_data[f'chiN{chiN_value}'] = data[:, i+1]

# Step 3: Assign chiN arrays to individual variables
gr1 = chiN_data[f'chiN{chiN[0]}']
gr2 = chiN_data[f'chiN{chiN[1]}']
gr3 = chiN_data[f'chiN{chiN[2]}']
gr4 = chiN_data[f'chiN{chiN[3]}']

# Clean gr
gr_clean1 = gr1 - hr1
gr_clean2 = gr1 - hr2
gr_clean3 = gr1 - hr3

#plotting
plt.plot(r, gr1, '-', label=r'$\chi$N=%d'%chiN[0], markersize=2, color=cmap1(0.23+float(0)/5))
# plt.plot(r, gr2, '-', label=r'$\chi$N=%d'%chiN[1], markersize=2, color=cmap1(0.23+float(1)/5))
plt.plot(r, gr_clean1,'--', label=r'$\chi$N=%d skew-normal fit'%chiN[1], markersize=2, color=cmap1(0.23+float(1)/5))
plt.plot(r, gr3, '-', label=r'$\chi$N=%d'%chiN[2], markersize=2, color=cmap1(0.23+float(2)/5))
plt.plot(r, gr_clean2, '--', label=r'$\chi$N=%d skew-normal fit'%chiN[2], markersize=2, color=cmap1(0.23+float(2)/5))
plt.plot(r, gr4, '-', label=r'$\chi$N=%d'%chiN[3], markersize=2, color=cmap1(0.23+float(3)/5))
plt.plot(r, gr_clean3, '--', label=r'$\chi$N=%d skew-normal fit'%chiN[3], markersize=2, color=cmap1(0.23+float(3)/5))
# plt.xlim(0, 20)
# plt.ylim(0,0.2)
plt.xlabel('r[b]')
plt.ylabel('g(r)')
plt.legend()
plt.savefig('gr_compare_fit.pdf')
plt.show()

# Save the cleaned gr data
np.savetxt("RDFAB_chiN%d.txt"%chiN[0], np.column_stack((r, gr1)), fmt='%.5f')
np.savetxt("RDFAB_chiN%d.txt"%chiN[1], np.column_stack((r, gr_clean1)), fmt='%.5f')
np.savetxt("RDFAB_chiN%d.txt"%chiN[2], np.column_stack((r, gr_clean2)), fmt='%.5f')
np.savetxt("RDFAB_chiN%d.txt"%chiN[3], np.column_stack((r, gr_clean3)), fmt='%.5f')
