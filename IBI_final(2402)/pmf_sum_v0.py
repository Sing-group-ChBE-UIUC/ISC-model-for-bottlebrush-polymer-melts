import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
import re
font = {'family': 'serif', 'serif': ['Computer Modern'], 'size': 16}
plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.rcParams['text.usetex'] = False

r1, pmf1, mf1 = np.loadtxt(sys.argv[1], unpack=True)
r2, pmf2, mf2 = np.loadtxt(sys.argv[2], unpack=True)
r3, pmf3, mf3 = np.loadtxt(sys.argv[3], unpack=True)
r4, pmf4, mf4 = np.loadtxt(sys.argv[4], unpack=True)
index, r5, pmf5, mf5 = np.loadtxt(sys.argv[5], unpack=True)
# Extracting chiN value from the file name
chiN = []

for filename in sys.argv[1:]:
    # Step 1: Extract \chi N value from the filename using regex
    match = re.search(r'chiN(\d+)', filename)
    if match:
        chiN_value = match.group(1)
        chiN.append(chiN_value)
    else:
        print(f"Filename {filename} does not match expected pattern.")
        continue


# Visualization

r=[r1,r2,r3, r4]#,r5,r6]
pmf=[pmf1,pmf2,pmf3,pmf4]#,pmf5,pmf6]
mf=[mf1,mf2,mf3,mf4]#,mf5,mf6]
label2=chiN

cmap1 = plt.get_cmap('Blues')
cmap2 = plt.get_cmap('Oranges')

fig, ax1 = plt.subplots()
#label=[r'$\chi$N=1', r'$\chi$N=13', r'$\chi$N=20', r'$\chi$N=20', r'$\chi$N=30']
colormap = ['red','blue','green','orange','purple',]
for i in range(4):
   ax1.plot(r[i],pmf[i], label=r'$\chi$N=%s'%label2[i], color = cmap1(0.3+(float(i))/4))
ax1.plot(r5, pmf5, label=r'$\chi$N=0', color='black', linestyle='--')


# plt.plot
plt.xlabel('r[b]')
plt.ylabel(r'V(r)/$k_{B}$T')
plt.ylim(-0.5,4.5)
plt.xlim(0, 15)
plt.legend()
plt.savefig('pot_chiN.pdf')
plt.show()
plt.close()
nbin = 440
r1 = r1[0:nbin]
pmf5 = pmf5[0:nbin]
pmf1 = pmf1[0:nbin]
pmf2 = pmf2[0:nbin]
pmf3 = pmf3[0:nbin]
pmf4 = pmf4[0:nbin]

header = f"r1ab chiN0 chiN{chiN[0]} chiN{chiN[1]} chiN{chiN[2]} chiN{chiN[3]}"
with open("pot_final.txt", "w") as f:
    f.write(header)
    f.write("\n")
    np.savetxt(f, np.column_stack((r1, pmf5, pmf1, pmf2, pmf3, pmf4)), fmt='%.5f')
    f.close()


dev1 = np.abs(pmf1 - pmf1)
dev2 = np.abs(pmf2 - pmf1)
dev3 = np.abs(pmf3 - pmf1)
dev4 = np.abs(pmf4 - pmf1)


# # # save data
header = f"r1ab chiN{chiN[1]} chiN{chiN[2]} chiN{chiN[3]}"
with open("devpot_final.txt", "w") as f:
    f.write(header)
    f.write("\n")
    np.savetxt(f, np.column_stack((r1, dev2, dev3, dev4)), fmt='%.5f')
    f.close()
# # ax2 = ax1.twinx()
fig, ax1 = plt.subplots()
ax1.plot(r1, dev1, label=r'$\chi$N=%s'%label2[0], color='red', linewidth=2)
ax1.plot(r1, dev2, label=r'$\chi$N=%s'%label2[1], color='blue', linewidth=2)
ax1.plot(r1, dev3, label=r'$\chi$N=%s'%label2[2], color='green', linewidth=2)
ax1.plot(r1, dev4, label=r'$\chi$N=%s'%label2[3], color='orange', linewidth=2)



plt.xlabel('r[b]')
# ax1.set_ylabel(r'$V(r)/k_{B}T$')
ax1.set_ylabel(r'$h_V(r)$')
# plt.xlim(1, 13)
# plt.ylim(0, 1)
plt.legend()
plt.savefig('dev.pdf')
plt.show()
plt.close()