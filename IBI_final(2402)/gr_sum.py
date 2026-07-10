import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import re

font = {'family': 'serif', 'serif': ['Computer Modern'], 'size': 15}
plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.rcParams['text.usetex'] = False
plt.figure(figsize=(9,5))
ax1=plt.subplot(1,1,1)

r1ab, grab1 = np.loadtxt(sys.argv[1], unpack=True)
r2ab, grab2 = np.loadtxt(sys.argv[2], unpack=True)
r3ab, grab3 = np.loadtxt(sys.argv[3], unpack=True)
r4ab, grab4 = np.loadtxt(sys.argv[4], unpack=True)
r1aa, graa1 = np.loadtxt(sys.argv[5], unpack=True)
# r1ab, grab1, mf = np.loadtxt(sys.argv[1], unpack=True)
# r2ab, grab2, mf = np.loadtxt(sys.argv[2], unpack=True)
# r3ab, grab3, mf = np.loadtxt(sys.argv[3], unpack=True)
# r4ab, grab4, mf = np.loadtxt(sys.argv[4], unpack=True)
# r1aa, graa1, mf = np.loadtxt(sys.argv[5], unpack=True)

r2aa, graa2 = np.loadtxt(sys.argv[6], unpack=True)
r3aa, graa3 = np.loadtxt(sys.argv[7], unpack=True)
r4aa, graa4 = np.loadtxt(sys.argv[8], unpack=True)
ro, gro = np.loadtxt(sys.argv[9], unpack=True)
D = float(sys.argv[10])
# d1 = float(sys.argv[6])
# d2 = float(sys.argv[7])
# d3 = float(sys.argv[8])
# d4 = float(sys.argv[9])
# d5 = float(sys.argv[10])


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

cmap1 = plt.get_cmap('Blues')
cmap2 = plt.get_cmap('Oranges')
cmap3 = ['darkviolet','cornflowerblue','mediumspringgreen', 'orange', 'orangered']
# raa=[r1aa/D,r2aa/D,r3aa/D,r4aa/D]
rab=[r1ab/D,r2ab/D,r3ab/D,r4ab/D]
grab=[grab1,grab2,grab3,grab4]
graa=[graa1,graa2,graa3,graa4]
# colormap1 = ['black','red','blue','green','orange','purple']

ax1.plot(ro/D, gro, label=r'$\chi$N=0',color='black', linestyle='-', linewidth=2.5)

for i in range(4):
    ax1.plot(rab[i],graa[i],label=r'$AA\ \chi$N=%s'%chiN[i], linestyle='-', color = cmap1(0.23+float(4-i)/4))
    
for i in range(4):
    ax1.plot(rab[i],grab[i], label=r'$AB\ \chi$N=%s'%chiN[i], linestyle='-', color = cmap2(0.23+float(4-i)/4))
  
ax1.set_xlim(0, 2.0)
ax1.set_ylim(0, 1.5)
ax1.xaxis.set_tick_params(pad=10)
ax1.yaxis.set_tick_params(pad=5)
plt.xlabel('r/D')
plt.ylabel('g(r)')
plt.legend(loc='lower right', fontsize=11, ncol=2)
plt.savefig('gr_figure.pdf', bbox_inches='tight')
plt.show()
plt.close()

# header = f"r1ab chiN0 chiN{chiN[0]} chiN{chiN[1]} chiN{chiN[2]} chiN{chiN[3]}"
# with open("gr_sum_final.txt", "w") as f:
#     f.write(header)
#     f.write("\n")
#     np.savetxt(f, np.column_stack((r1ab, gro, grab1, grab2, grab3, grab4)), fmt='%.5f')
#     f.close()

# ## homo compare
# f1_r=[r1ab, r2ab, r3ab, r4ab, r1aa]
# f1_gr=[grab1, grab2, grab3, grab4, graa1]
# D=[d1, d2, d3, d4, d5]
# # f2_r=[r4ab, r1aa, r2aa]
# # f2_gr=[grab4, graa1, graa2]
# # f5_r=[r3aa, r4aa, ro]
# # f5_gr=[graa3, graa4, gro]
# Nsc=[2, 4, 8, 14, 20]

# for i in range(5):
#     # ax1.plot(f1_r[i],f1_gr[i],label=r'$N_{sc}$=%d'%Nsc[i], linestyle='-', color = cmap3[i], linewidth=2.5)
#     ax1.plot(f1_r[i]/D[i],f1_gr[i],label=r'$N_{sc}$=%d'%Nsc[i], linestyle='-', color = cmap3[i], linewidth=2.5)
# # ax1.plot(f1_r[4]/(10.8/11.8),f1_gr[4],label=r'$N_{sc}$=%d'%Nsc[4], linestyle='-', color = cmap3[4], linewidth=2.5)

    

# # ax1.set_xlim(0, 13)
# # ax1.set_xlim(0, 2)
# # ax1.set_ylim(0, 1.5)
# ax1.xaxis.set_tick_params(pad=10)
# ax1.yaxis.set_tick_params(pad=5)
# plt.xlabel('r/D')
# plt.ylabel('g(r)')
# plt.legend(loc='lower right', fontsize=15, ncol=1)
# plt.savefig('gr_homofigure_f5.pdf', bbox_inches='tight')
# plt.show()
# plt.close()


# for i in range(3):
#     ax1.plot(f2_r[i],f2_gr[i],label=r'$N_{sc}$=%d'%Nsc[i], linestyle='-', color = cmap3(0.23+float(4-i)/4))
    

# ax1.set_xlim(0, 30)
# ax1.set_ylim(0, 1.5)
# ax1.xaxis.set_tick_params(pad=10)
# ax1.yaxis.set_tick_params(pad=5)
# plt.xlabel('r[b]')
# plt.ylabel('g(r)')
# plt.legend(loc='lower right', fontsize=11, ncol=2)
# plt.savefig('gr_homofigure_f2.pdf', bbox_inches='tight')
# plt.show()


# for i in range(3):
#     ax1.plot(f5_r[i],f5_gr[i],label=r'$N_{sc}$=%d'%Nsc[i], linestyle='-', color = cmap3(0.23+float(4-i)/4))
    

# ax1.set_xlim(0, 30)
# ax1.set_ylim(0, 1.5)
# ax1.xaxis.set_tick_params(pad=10)
# ax1.yaxis.set_tick_params(pad=5)
# plt.xlabel('r[b]')
# plt.ylabel('g(r)')
# plt.legend(loc='lower right', fontsize=11, ncol=2)
# plt.savefig('gr_homofigure_f5.pdf', bbox_inches='tight')
# plt.show()
# plt.close()