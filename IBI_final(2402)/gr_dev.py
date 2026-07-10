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
# r5ab, grab5 = np.loadtxt(sys.argv[5], unpack=True)
# r6ab, grab6 = np.loadtxt(sys.argv[6], unpack=True)
# r4aa, graa4 = np.loadtxt(sys.argv[7], unpack=True)
# r4ab, grab4 = np.loadtxt(sys.argv[8], unpack=True)
# r5, graa5 = np.loadtxt(sys.argv[9], unpack=True)
# r5, grab5 = np.loadtxt(sys.argv[10], unpack=True)
# r6, graa6 = np.loadtxt(sys.argv[11], unpack=True)
# r6, grab6 = np.loadtxt(sys.argv[12], unpack=True)
cmap1 = plt.get_cmap('Blues')
cmap2 = plt.get_cmap('Oranges')
# raa=[r1aa,r2aa,r3aa,r4aa]#,r5]#,r6]
rab=[r1ab,r2ab,r3ab,r4ab]#,r5ab,r6ab]

# Extracting chiN value from the file name
chiN = []
for i in range(1,4):
    for filename in sys.argv[1:]:
        # Step 1: Extract \chi N value from the filename using regex
        match = re.search(r'chiN(\d+)', filename)
        if match:
            chiN_value = match.group(1)
            chiN.append(chiN_value)
        else:
            print(f"Filename {filename} does not match expected pattern.")
            continue
        
# graa=[graa1,graa2,graa3,graa4]#,graa5]#,graa6]
grab=[grab1,grab2,grab3,grab4]#, grab5,grab6]
colormap1 = ['black','red','blue','green','orange','purple']
dev1 = grab1-grab2
dev2 = grab1-grab3
dev3 = grab1-grab4
#dev4 = grab1-grab5
#dev5 = grab1-grab6
gr_dev = [dev1, dev2, dev3]#, dev4, dev5]

# ro, gro = np.loadtxt(sys.argv[9], unpack=True)
# plt.plot(ro, gro, label=r'$\chi$N=0',color='black', linestyle='-')

# for i in range(4):
#     plt.plot(raa[i],graa[i],label=r'$g_{AA}\ \chi$N=%d'%chiN[i], linestyle='-', color = cmap1(0.23+float(i)/4))
    
for i in range(3):
    plt.plot(rab[i],gr_dev[i], label=r'$\chi$N=%s'%chiN[i+1], linestyle='-', color = cmap2(0.23+float(i)/4))

header = f"r1ab chiN{chiN[1]} chiN{chiN[2]} chiN{chiN[3]}"
with open("dev_org_final.txt", "w") as f:
    f.write(header)
    f.write("\n")
    np.savetxt(f, np.column_stack((r1ab, gr_dev[0], gr_dev[1], gr_dev[2])), fmt='%.5f')
    f.close()
plt.xlim(0, 15)
# plt.ylim(0, 0.16)
plt.xlabel('r[b]')
plt.ylabel('h(r)')
plt.legend(loc='upper right', fontsize=11, ncol=1)
plt.savefig('gr_dev.pdf')
plt.show()
plt.close()

