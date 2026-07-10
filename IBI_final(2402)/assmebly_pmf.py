import numpy as numpy
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


r, gr = numpy.loadtxt(sys.argv[1], unpack=True)

def power(x, a, b):
    return a*(x/5.4)**b

popt, pcov = curve_fit(power, r, gr)
print(popt)
print(pcov)
new_r = numpy.linspace(r[0], r[-1], 100)
new_gr = power(new_r, *popt)
plt.plot(r, gr, label='g(r)')
plt.plot(new_r, new_gr, label='fit')
plt.legend()
plt.show()
plt.close()