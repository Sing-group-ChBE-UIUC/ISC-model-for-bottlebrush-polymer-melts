import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import skewnorm

font = {'family': 'serif', 'serif': ['Computer Modern'], 'size': 16}
plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.rcParams['text.usetex'] = False
plt.figure(figsize=(8,6))
ax1=plt.subplot(1,1,1)

# Initialize lists and arrays to store data
nsc = [2,4,8,14,20]
f = [1,2,5]
for side in range(len(nsc)):
    chiN = []
    r1ab = []
    chiN_data = {}
    skew_fit_data = {}
    # print(side, 3*side, 3*side+1, 3*side+2)
    with open(sys.argv[side+1], "r") as file:
        header = file.readline().strip().split()
        chiN = [int(col.replace('chiN', '')) for col in header[1:]]
        data = np.loadtxt(file)
        r = data[:, 0]
        for i, chiN_value in enumerate(chiN):
            chiN_data[f'chiN{chiN_value}'] = data[:, i+1]

    # ro, gro = np.loadtxt(sys.argv[3*side+2], unpack=True)
    # r1, gr1 = np.loadtxt(sys.argv[3*side+3], unpack=True)
    # with open(sys.argv[2], "r") as file:
    #     header = file.readline().strip().split()
    #     file.readline().strip().split()
    #     file.readline().strip().split()
    #     file.readline().strip().split()
    #     data = np.loadtxt(file)
    #     r = data[:, 0]
    #     for i, chiN_value in enumerate(chiN):
    #         skew_fit_data[f'chiN{chiN_value}'] = data[:, i+1]

    # D=float(sys.argv[2])
    # Assign chiN arrays to individual variables
    hr1 = chiN_data[f'chiN{chiN[0]}']
    hr2 = chiN_data[f'chiN{chiN[1]}']
    hr3 = chiN_data[f'chiN{chiN[2]}']
    hr4= chiN_data[f'chiN{chiN[3]}']
    # # Assign skew_fit arrays to individual variables
    # fit1 = skew_fit_data[f'chiN{chiN[0]}']
    # fit2 = skew_fit_data[f'chiN{chiN[1]}']
    # fit3 = skew_fit_data[f'chiN{chiN[2]}']


    cmap1 = plt.get_cmap('CMRmap')

    h_r=[hr1, hr2, hr3, hr4]
    # fit=[fit1, fit2, fit3]


    # for i in range(3):

    #     # plt.plot(r/D, h_r[i], 'o', label=r'$\chi$N=%d'%chiN[i], alpha=0.5, markersize=5, markerfacecolor = 'none', markeredgecolor=cmap1(0.23+float(i)/4), markeredgewidth=1.5)
    #     plt.plot(r/D, h_r[i], '-', label=r'$\chi$N=%d'%chiN[i], alpha=0.5,color=cmap1(0.23+float(i)/4), linewidth=3.5)
    # for i in range(3):
    #     plt.plot(r/D, fit[i], '--', label=r'$\chi$N=%d skew-normal '%chiN[i], color=cmap1(0.23+float(i)/4), linewidth=2.5)

        
    # plt.xlim(0, 2)
    # plt.ylim(0,0.3)
    # plt.xlabel('r/D')
    # plt.ylabel('h(r)')
    # plt.legend(fontsize=15)
    # plt.savefig('fit_figure.pdf', bbox_inches='tight')
    # plt.show()
    # plt.close()

    # ## This section is to overlap the data
    cmap3 = plt.get_cmap('CMRmap')
    # len1 = len(ro)
    # len2 = len(r1)
    # if(len1>len2):
    #     r = r[:len2]
    #     ro = ro[:len2]
    #     gro = gro[:len2]
    #     hr1 = hr1[:len2]
    #     hr2 = hr2[:len2]
    #     hr3 = hr3[:len2]
        
    # else: 
    #     if (len1<len2):
    #         r = r[:len1]
    #         r1 = r1[:len1]
    #         gr1 = gr1[:len1]
    #         hr1 = hr1[:len1]
    #         hr2 = hr2[:len1]
    #         hr3 = hr3[:len1]
    # dev0= gro-gr1
    # dev1 = hr1+dev0
    # dev2 = hr2+dev0
    # dev3 = hr3+dev0
    dev0 = hr1
    dev1 = hr2
    dev2 = hr3
    dev3 = hr4


    A_hight0=max(dev0[50:])
    A_hight1=max(dev1[50:])
    A_hight2=max(dev2[30:])
    A_hight3=max(dev3[30:])
    xi0=r[np.argmax(dev0[50:])]
    xi1=r[np.argmax(dev1[50:])]
    xi2=r[np.argmax(dev2[30:])]
    xi3=r[np.argmax(dev3[30:])]
    m1 = 0
    m2 = 0
    m3 = 0
    m4 = 0
    n1 = 0
    n2 = 0
    n3 = 0
    n4 = 0
    if side == 0:
        m3 = 0.6
        m4 = 0.6
    if side == 1:
        m3 = 0.5
        m4 = 0.6
        n3 = 0.
        n4 = 0.
        
    if side == 2:
        m3 = 0.4
        m4 -0.3
    if side == 3:
        m3 = -0.3
        m4 = 0.3

        
    if side == 4:
        m3 = 0.2
        m4 = 0.2
   
      

    # ax1.plot(r/(xi0+m1), dev0/(A_hight0+n1), label=r'$\chi$N=1', color=cmap3(0.23+float(0)/5), linewidth=2)
    # ax1.plot(r/(xi1+m2), dev1/(A_hight1+n2), label=r'$\chi$N=%s'%chiN[0], color=cmap3(0.23+float(1)/5), linewidth=2)
    # ax1.plot(r/(xi2+m3), dev2/(A_hight2+n3), label=r'$\chi$N=%s'%chiN[1], color=cmap3(0.23+float(2)/5), linewidth=2)
    # ax1.plot(r/(xi3+m4), dev3/(A_hight3+n4), label=r'$\chi$N=%s'%chiN[2], color=cmap3(0.23+float(3)/5), linewidth=2)
    ax1.plot(r/(xi2+m3), dev2/(A_hight2+n3), label=r'$N_{sc}$=%s'%nsc[side], color=cmap3(float(side)/5), linewidth=2)
    # ax1.plot(r/(xi2+m3), dev2/(A_hight2+n3), label=r'$f$=%s'%f[side], color=cmap3(0.1+float(side)/3), linewidth=2)
    ax1.plot(r/(xi3+m4), dev3/(A_hight3+n4), alpha = 0.5, color=cmap3(float(side)/5), linewidth=2)
    
    # if side == 4:
    #     ax1.plot(r/(xi2+m3), dev2/(A_hight2+n3), alpha=0.5, color=cmap3(float(side)/5), linewidth=2)
    #     ax1.plot(r/(xi3+m4), dev3/(A_hight3+n4), label=r'$N_{sc}$=%s'%nsc[side], color=cmap3(float(side)/5), linewidth=2)
    # else:
    #     ax1.plot(r/(xi2+m3), dev2/(A_hight2+n3), label=r'$N_{sc}$=%s'%nsc[side], color=cmap3(float(side)/5), linewidth=2)
    #     ax1.plot(r/(xi3+m4), dev3/(A_hight3+n4), alpha = 0.5, color=cmap3(float(side)/5), linewidth=2)
        
    # ax1.plot(r/(xi2+0.1), dev3/A_hight3, alpha = 0.5, color=cmap3(0.23+float(0)/5), linewidth=2)
    print("Nsc = %d"%(nsc[side]), f"{xi0+m1:.3f}", f"{xi1+m2:.3f}", f"{xi2+m3:.3f}", f"{xi3+m4:.3f}")
    print("Nsc = %d"%(nsc[side]), f"{A_hight0 + n1:.3f}", f"{A_hight1 + n2:.3f}", f"{A_hight2 + n3:.3f}", f"{A_hight3 + n4:.3f}")


plt.xlim(0,4)
plt.ylim(0,2)
plt.xlabel(r'$r/\xi$')
plt.ylabel(r'$h_{V}(r)/A$')
plt.legend(fontsize=15)
plt.savefig('potdev_overlap_f5.pdf', bbox_inches='tight')
plt.show()
plt.close()
# print(f"{xi0+m1:.3f}", f"{xi1+m2:.3f}", f"{xi2+m3:.3f}", f"{xi3+m4:.3f}")
# print(xi0+m1, xi1+m2, xi2+m3, xi3+m4)
# print(A_hight0+n1, A_hight1+n2, A_hight2+n3, A_hight3+n4)
# print(f"{A_hight0 + n1:.3f}", f"{A_hight1 + n2:.3f}", f"{A_hight2 + n3:.3f}", f"{A_hight3 + n4:.3f}")
