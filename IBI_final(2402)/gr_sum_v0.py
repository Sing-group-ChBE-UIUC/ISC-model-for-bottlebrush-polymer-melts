import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import re

font = {'family': 'serif', 'serif': ['Computer Modern'], 'size': 16}
plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.rcParams['text.usetex'] = False

r1ab, grab1 = np.loadtxt(sys.argv[1], unpack=True)
r2ab, grab2 = np.loadtxt(sys.argv[2], unpack=True)
r3ab, grab3 = np.loadtxt(sys.argv[3], unpack=True)
r4ab, grab4 = np.loadtxt(sys.argv[4], unpack=True)
# r3aa, graa3 = np.loadtxt(sys.argv[5], unpack=True)
# r3ab, grab3 = np.loadtxt(sys.argv[6], unpack=True)
# r4aa, graa4 = np.loadtxt(sys.argv[7], unpack=True)
# r4ab, grab4 = np.loadtxt(sys.argv[8], unpack=True)
# r5, graa5 = np.loadtxt(sys.argv[9], unpack=True)
# r5, grab5 = np.loadtxt(sys.argv[10], unpack=True)
# r6, graa6 = np.loadtxt(sys.argv[11], unpack=True)
# r6, grab6 = np.loadtxt(sys.argv[12], unpack=True)

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
# raa=[r1aa,r2aa,r3aa,r4aa]#,r5]#,r6]
rab=[r1ab,r2ab,r3ab,r4ab]#,r5]#,r6]
# chiN=[1, 17, 33, 48]
# graa=[graa1,graa2,graa3,graa4]#,graa5]#,graa6]
grab=[grab1,grab2,grab3,grab4]#,grab5]#,grab6]
colormap1 = ['black','red','blue','green','orange','purple']
print(chiN[1])

# ro, gro = np.loadtxt(sys.argv[9], unpack=True)
# plt.plot(ro, gro, label=r'$\chi$N=0',color='black', linestyle='-')

# for i in range(4):
#     plt.plot(raa[i],graa[i],label=r'$g_{AA}\ \chi$N=%d'%chiN[i], linestyle='-', color = cmap1(0.23+float(i)/4))
    
for i in range(4):
    plt.plot(rab[i],grab[i], label=r'$\chi$N=%s'%chiN[i], linestyle='-', color = cmap2(0.23+float(i)/4))
  
# r, gro = np.loadtxt(sys.argv[9], unpack=True)
# plt.plot(r, gro, label=r'$AA\chi$N=0',color='black', linestyle='-')
# plt.plot(r[5],graa[5],label=r'$\chi$N=20', color=colormap[5], linestyle='--')
# plt.plot(r[5],grab[5], color=colormap[5])
# plt.xlim(0, 20)
# plt.ylim(0, 1.6)
plt.xlabel('r[b]')
plt.ylabel('g(r)')
plt.legend(loc='lower right', fontsize=11, ncol=2)
plt.savefig('gr.pdf')
plt.show()
plt.close()
header = f"r1ab chiN{chiN[0]} chiN{chiN[1]} chiN{chiN[2]} chiN{chiN[3]}"
with open("gr_org_final.txt", "w") as f:
    f.write(header)
    f.write("\n")
    np.savetxt(f, np.column_stack((r1ab, grab1, grab2, grab3, grab4)), fmt='%.5f')
    f.close()
# np.savetxt("gr_org_final.txt", np.column_stack((r1ab, grab1, grab2, grab3, grab4)), fmt='%.5f', header=header)

