import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
import re
font = {'family': 'serif', 'serif': ['Computer Modern'], 'size': 16}
plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.rcParams['text.usetex'] = False
plt.figure(figsize=(8,5))
ax1=plt.subplot(1,1,1)

# # chiN0 potential by architecture
# r0, pmf0= np.loadtxt(sys.argv[1], unpack=True)
# r1, pmf1= np.loadtxt(sys.argv[2], unpack=True)
# r2, pmf2= np.loadtxt(sys.argv[3], unpack=True)
# D1 = float(sys.argv[4])
# D2 = float(sys.argv[5])
# D3 = float(sys.argv[6])
# cmap3 = ['darkviolet','cornflowerblue','mediumspringgreen', 'orange', 'orangered']

# ax1.plot(r0/D1, pmf0, label=r'$N_{sc}$=2', color=cmap3[0], linestyle='-', linewidth=1.5)
# ax1.plot(r1/D2, pmf1, label=r'$N_{sc}$=4', color=cmap3[1], linestyle='-', linewidth=1.5)
# ax1.plot(r2/D3, pmf2, label=r'$N_{sc}$=8', color=cmap3[2], linestyle='-', linewidth=1.5)
# plt.xlim(0,1.5)
# plt.ylim(-1,4)
# plt.xlabel('r/D')
# plt.ylabel(r'V(r)/$k_{B}$T')
# plt.legend()
# plt.savefig('V(r)_f5.pdf', bbox_inches='tight')
# plt.show()
# plt.close()


# Get pot profile by chiN based on pot0.table

# index, r0, pmf0, mf0 = np.loadtxt(sys.argv[1], unpack=True)
r0, pmf0= np.loadtxt(sys.argv[1], unpack=True)
D=float(sys.argv[3])
# Extracting chiN value from the file name
chiN = []
chiN_data = {}

with open(sys.argv[2], "r") as file:
    header = file.readline().strip().split()
    chiN = [int(col.replace('chiN', '')) for col in header[1:]]
    data = np.loadtxt(file)
    r = data[:, 0]
    for i, chiN_value in enumerate(chiN):
        chiN_data[f'chiN{chiN_value}'] = data[:, i+1]


# Visualization

# pmf0 = pmf0[:len(r)]
dev0 = chiN_data[f'chiN{chiN[0]}']
dev1 = chiN_data[f'chiN{chiN[1]}']
dev2 = chiN_data[f'chiN{chiN[2]}']
dev3 = chiN_data[f'chiN{chiN[3]}']
# dev0 = dev0[:len(pmf0)]
# dev1 = dev1[:len(pmf0)]
# dev2 = dev2[:len(pmf0)]
# dev3 = dev3[:len(pmf0)]
# r= r[:len(pmf0)]
dev=[dev1,dev2,dev3]#,dev4]#,dev5,dev6]

# cmap1 = plt.get_cmap('Blues')

# pmf1 = pmf0+dev0
# pmf2 = pmf0+dev1
# pmf3 = pmf0+dev2
# pmf4 = pmf0+dev3

# pmf=[pmf1,pmf2,pmf3,pmf4]#,pmf5,pmf6]

A_hight0=max(dev0[20:])
A_hight1=max(dev1[20:])
A_hight2=max(dev2[20:])
A_hight3=max(dev3[20:])
xi0=r[np.argmax(dev0[20:])]
xi1=r[np.argmax(dev1[20:])]
xi2=r[np.argmax(dev2[20:])]
xi3=r[np.argmax(dev3[20:])]
# scale = 5
# pmf0=pmf0[:len(dev2)]
# r0=r0[:len(dev2)]
# r0_norm = r0/D
# new_pot = pmf0+dev2*scale

# # bring section for extrapolation
# new_pot_predict = np.zeros(len(r0))
# indice = np.where((r0_norm>=0.5) & (r0_norm<=0.65))[0]
# r_test = r0_norm[indice[0]:indice[-1]]
# pot_test = new_pot[indice[0]:indice[-1]]

# # Find the best fit of curve
# def ext(x, a, b, c):
#     return a * np.exp(-b * x) + c
# # def log(x, a, b, c):
# #     return a * np.log(b * x) + c
# popt, pcov = curve_fit(ext, r_test, pot_test)
# f1 = ext(r0_norm, *popt)

# ax1.plot(r0_norm[:indice[-1]-3], f1[:indice[-1]-3], color="black", linewidth=2.)
# ax1.plot(r0_norm[indice[-1]-1:], new_pot[indice[-1]-1:], label=r'$\chi N_{\mathrm{LAM}}$', color="black", linewidth=2.)
# ax1.plot(r0_norm[:], new_pot[:], label=r'$\chi N_{\mathrm{LAM}}$', color="black", linewidth=2.)

for i in range(4):
   ax1.plot(r/D,pmf[3-i], label=r'$\chi$N=%s'%chiN[3-i], color = cmap1(0.3+(float(3-i))/4), linewidth=1.5)
ax1.plot(r/D, pmf0, label=r'$\chi$N=0', color='Orangered', linestyle='-', linewidth=1.5)

# plt.plot
plt.xlabel(r'r/D')
plt.ylabel(r'V(r)/$k_{B}$T')
plt.ylim(-1,4)
plt.xlim(0, 1.5)
plt.legend()
plt.savefig('pot_chiN.pdf')
plt.show()
plt.close()
































































# ## This section is to get potential difference

# # Extracting chiN value from the file name
# chiN = []
# chiN_data = {}

# with open(sys.argv[1], "r") as file:
#     header = file.readline().strip().split()
#     chiN = [int(col.replace('chiN', '')) for col in header[1:]]
#     data = np.loadtxt(file)
#     r = data[:, 0]
#     for i, chiN_value in enumerate(chiN):
#         chiN_data[f'chiN{chiN_value}'] = data[:, i+1]

# cmap3 = plt.get_cmap('CMRmap')

# D=float(sys.argv[2])
# # Visualization

# pmf0 = chiN_data[f'chiN{chiN[0]}']
# pmf1 = chiN_data[f'chiN{chiN[1]}']
# pmf2 = chiN_data[f'chiN{chiN[2]}']
# pmf3 = chiN_data[f'chiN{chiN[3]}']
# pmf4 = chiN_data[f'chiN{chiN[4]}']
# pmf=[pmf1,pmf2,pmf3,pmf4]#,pmf5,pmf6]


# dev1 = np.abs(pmf1 - pmf0)
# dev2 = np.abs(pmf2 - pmf0)
# dev3 = np.abs(pmf3 - pmf0)
# dev4 = np.abs(pmf4 - pmf0)


# # # # # save data
# header = f"r1ab chiN{chiN[1]} chiN{chiN[2]} chiN{chiN[3]} chiN{chiN[4]}"
# with open("devpot_final.txt", "w") as f:
#     f.write(header)
#     f.write("\n")
#     np.savetxt(f, np.column_stack((r, dev1, dev2, dev3, dev4)), fmt='%.5f')
#     f.close()
# # # # ax2 = ax1.twinx()
# r=r[3:len(dev1)]
# dev1=dev1[3:]
# dev2=dev2[3:]
# dev3=dev3[3:]
# dev4=dev4[3:]

# ax1.plot(r/D, dev1, label=r'$\chi$N=%s'%chiN[1], color=cmap3(0.23+float(0)/5), linewidth=2)
# ax1.plot(r/D, dev2, label=r'$\chi$N=%s'%chiN[2], color=cmap3(0.23+float(1)/5), linewidth=2)
# ax1.plot(r/D, dev3, label=r'$\chi$N=%s'%chiN[3], color=cmap3(0.23+float(2)/5), linewidth=2)
# ax1.plot(r/D, dev4, label=r'$\chi$N=%s'%chiN[4], color=cmap3(0.23+float(3)/5), linewidth=2)



# plt.xlabel('r/D')
# # ax1.set_ylabel(r'$V(r)/k_{B}T$')
# ax1.set_ylabel(r'$h_v(r)$')
# plt.xlim(0, 1.5)
# plt.ylim(0, 1.5)
# plt.legend()
# plt.savefig('dev_fig.pdf', bbox_inches='tight')
# plt.show()
# plt.close()