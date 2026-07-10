import numpy as np
import matplotlib.pyplot as plt
#from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
import sys

# r, pmf, mf = np.loadtxt(sys.argv[1], unpack = True)
# spline = UnivariateSpline(range(len(pmf)), pmf, s=0.1)
# pmf_sm = spline(range(len(pmf)))
# d=r[1]-r[0]
# mf=np.zeros(len(pmf_sm))
# for i in range(1, len(pmf_sm)-1):
#     mf[i] = -(pmf_sm[i+1]-pmf_sm[i-1])/(2*d)
# mf[0]=-(pmf_sm[1]-pmf_sm[0])/(d)
# mf[-1]=0    
# plt.plot(r, pmf_sm, label = 'smoothed')
# plt.ylim(-0.5,6)
# plt.show()
# filename="pot_sm.table"
# with open (filename, "w") as f:
#     for i in range(len(r)):
#         f.write("%d %lf %lf %lf\n"%(i+1, r[i], pmf_sm[i], mf[i]))
# f.close()



ro, gro= np.loadtxt(sys.argv[1], unpack = True)
#r, gr = np.loadtxt(sys.argv[2], unpack=True)
# index, ro, gro, mf = np.loadtxt(sys.argv[1], unpack=True)
spline = UnivariateSpline(range(len(gro)), gro, s=float(sys.argv[2]))
gr_sm = spline(range(len(gro)))
#gr_new = np.exp(-pmf)
#d=rp[1]-rp[0]
plt.plot(ro, gro, label = 'original')
#plt.plot(ro, gro, label='orginal')
plt.plot(ro, gr_sm, label='sm')
plt.legend()
#plt.ylim(-1,6)
#plt.xlim(0, 10)
plt.show()
filename="potsm.txt"
with open (filename, "w") as f:
    for i in range(len(ro)):
        # f.write("%d %lf %lf %lf\n"%(index[i], ro[i], gr_sm[i], mf[i]))
        f.write("%lf %lf\n"%(ro[i], gr_sm[i]))
f.close()

