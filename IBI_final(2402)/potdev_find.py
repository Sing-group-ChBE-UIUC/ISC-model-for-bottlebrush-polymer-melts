import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import re
font = {'family': 'serif', 'serif': ['Computer Modern'], 'size': 16}
plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.rcParams['text.usetex'] = False
plt.figure(figsize=(8,6))
ax1=plt.subplot(1,1,1)

# # Get maximum value of deviation
# Extracting chiN value from the file name
chiN = []
chiN_data = {}

with open(sys.argv[1], "r") as file:
    header = file.readline().strip().split()
    chiN = [int(col.replace('chiN', '')) for col in header[1:]]
    data = np.loadtxt(file)
    r = data[:, 0]
    for i, chiN_value in enumerate(chiN):
        chiN_data[f'chiN{chiN_value}'] = data[:, i+1]


# Visualization
dev1 = chiN_data[f'chiN{chiN[0]}']
dev2 = chiN_data[f'chiN{chiN[1]}']
dev3 = chiN_data[f'chiN{chiN[2]}']
dev4 = chiN_data[f'chiN{chiN[3]}']

dev=[dev1,dev2,dev3]#,dev4]#,dev5,dev6]

A_hight1=max(dev1[30:])
A_hight2=max(dev2[30:])
A_hight3=max(dev3[30:])
A_hight4=max(dev4[30:])

xi1=r[np.argmax(dev1[30:])]
xi2=r[np.argmax(dev2[30:])]
xi3=r[np.argmax(dev3[30:])]
xi4=r[np.argmax(dev4[30:])]

cmap3 = plt.get_cmap('CMRmap')

# D=float(sys.argv[2])
# Visualization

# r=r[3:len(dev1)]
# dev1=dev1[3:]
# dev2=dev2[3:]
# dev3=dev3[3:]
# dev4=dev4[3:]
# scale = 3
# ax1.plot(r/D, dev1, label=r'$\chi$N=%s'%chiN[0], color=cmap3(0.23+float(0)/5), linewidth=2)
# ax1.plot(r/D, dev2, label=r'$\chi$N=%s'%chiN[1], color=cmap3(0.23+float(1)/5), linewidth=2)
# ax1.plot(r/D, dev3, label=r'$\chi$N=%s'%chiN[2], color=cmap3(0.23+float(2)/5), linewidth=2)
# ax1.plot(r/D, dev4, label=r'$\chi$N=%s'%chiN[3], color=cmap3(0.23+float(3)/5), linewidth=2)
# ax1.plot(r/D, dev3*scale, label = r'$\chi N_{\mathrm{LAM}}$', linestyle='--', color="black", linewidth=2.)
ax1.plot(r/(xi1), dev1/(A_hight1), label=r'$\chi$N=%s'%chiN[0], color=cmap3(0.23+float(0)/5), linewidth=2)
ax1.plot(r/(xi2), dev2/(A_hight2), label=r'$\chi$N=%s'%chiN[1], color=cmap3(0.23+float(1)/5), linewidth=2)
ax1.plot(r/(xi2), dev3/A_hight3, label=r'$\chi$N=%s'%chiN[2], color=cmap3(0.23+float(2)/5), linewidth=2)
ax1.plot(r/(xi2), dev4/A_hight4, label=r'$\chi$N=%s'%chiN[3], color=cmap3(0.23+float(3)/5), linewidth=2)



# plt.xlabel(r'r/D')
# # ax1.set_ylabel(r'$V(r)/k_{B}T$')
ax1.set_ylabel(r'$h_v(r)$')
plt.xlim(0, 1.4)
plt.ylim(0, 3)
plt.legend()
plt.savefig('dev_f1.pdf', bbox_inches='tight')
plt.show()
plt.close()

# ax1=plt.subplot(1,1,1)
# r0, pot0 = np.loadtxt(sys.argv[3], unpack=True)
# pot0=pot0[:len(dev3)]
# r0=r0[:len(dev3)]
# r0_norm = r0/D
# new_pot = pot0+dev3*scale

# # bring section for extrapolation
# new_pot_predict = np.zeros(len(r0))
# indice = np.where((r0_norm>=0.5) & (r0_norm<=0.65))[0]
# r_test = r0_norm[indice[0]:indice[-1]]
# pot_test = new_pot[indice[0]:indice[-1]]
# # ax1.plot(r_test, pot_test, label='test', color="green", linewidth=2.)

# # Find the best fit of curve
# def ext(x, a, b, c):
#     return a * np.exp(-b * x) + c
# # def log(x, a, b, c):
# #     return a * np.log(b * x) + c
# popt, pcov = curve_fit(ext, r_test, pot_test)
# f1 = ext(r0_norm, *popt)

# # new_pot_predict[indice[-1]:] = new_pot[indice[-1]:]
# ax1.plot(r0_norm[:indice[-1]], f1[:indice[-1]], label=r'exponential fit', color="orange", linewidth=2.)
# ax1.plot(r0/D, new_pot, linestyle='--', label=r'$\chi N_{\mathrm{LAM}}$', color="black", linewidth=2.)
# ax1.plot(r0/D, pot0, label=r'$\chi N=0$', color="Orangered", linewidth=2.)

# # ax1.plot(r0/D, new_pot_predict, label=r'corr\ $\chi N_{\mathrm{LAM}}$', color="black", linewidth=2.)
# plt.xlabel(r'r/D')
# # ax1.set_ylabel(r'$V(r)/k_{B}T$')
# ax1.set_ylabel(r'V(r)/$k_{B}$T')
# plt.xlim(0, 1.4)
# plt.ylim(-1, 5)
# plt.legend()
# plt.savefig('new_pot.pdf', bbox_inches='tight')
# plt.show()
# plt.close()

# print(chiN[0], chiN[1], chiN[2], chiN[3])
# print(A_hight1, A_hight2, A_hight3, A_hight4)
# print(xi1, xi2, xi3, xi4)