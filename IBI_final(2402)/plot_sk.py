#!/usr/bin/env python

import math
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as ticr
import pandas as pd
import sys
import matplotlib.pyplot as plt

x1, y1 = np.loadtxt(sys.argv[1],unpack=True)
x2, y2 = np.loadtxt(sys.argv[2],unpack=True)
x3, y3 = np.loadtxt(sys.argv[3],unpack=True)


num_bins = 1000
sum_y1 = np.histogram(x1,bins=num_bins,weights=y1)[0]
sum_y2 = np.histogram(x2,bins=num_bins,weights=y2)[0]
sum_y3 = np.histogram(x3,bins=num_bins,weights=y3)[0]
count_y1 = np.histogram(x1,bins=num_bins)[0]
count_y2 = np.histogram(x2,bins=num_bins)[0]
count_y3 = np.histogram(x3,bins=num_bins)[0]
xbin1 = np.histogram(x1,bins=num_bins)[1][1:][sum_y1 != 0.0]
xbin2 = np.histogram(x2,bins=num_bins)[1][1:][sum_y2 != 0.0]
xbin3 = np.histogram(x3,bins=num_bins)[1][1:][sum_y3 != 0.0]
ave_y1 = sum_y1[sum_y1 != 0.0]/count_y1[count_y1 != 0]*20
ave_y2 = sum_y2[sum_y2 != 0.0]/count_y2[count_y2 != 0]*20
ave_y3 = sum_y3[sum_y3 != 0.0]/count_y3[count_y3 != 0]*20

xlist=[xbin1, xbin2, xbin3]
ylist=[ave_y1, ave_y2, ave_y3]

markers = ['s']*14
cmap = mpl.cm.get_cmap('rainbow', 7)
# cmap = ['aquamarine']

# xdata1 = xbin1*0.01317
# ydata1 = ave_y1
# plt.plot(xdata1[2:], ydata1[2:], marker=markers[0], mec=cmap[0], ms=1.5, color=cmap[0])
plt.plot(xlist[0][2:], ylist[0][2:], marker=markers[0], mec=cmap(0), ms=1.5, color=cmap(0), label = r'$N_{sc}=2$')#, label = 'f2nsc20')
plt.plot(xlist[1][2:], ylist[1][2:], marker=markers[1], mec=cmap(1), ms=1.5, color=cmap(1), label = r'$N_{sc}=8$')#, label = 'f5nsc8-f2nsc20')
plt.plot(xlist[2][2:], ylist[2][2:], marker=markers[2], mec=cmap(2), ms=1.5, color=cmap(2), label = r'$N_{sc}=20$')#, label = 'f5nsc20-f2nsc20')

plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$S(q)\ [a.u.]$')
plt.xlabel(r'$q\ [b^{-1}]$')
# plt.xlim(0.03,10)
# plt.ylim(20, 22)
plt.legend()
plt.savefig('sk.png')
plt.show()
plt.close()

