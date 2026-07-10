import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import skewnorm

r, hr1, hr2, hr3= np.loadtxt(sys.argv[1], unpack=True)
# r, hr2= np.loadtxt(sys.argv[2], unpack=True)

# # Initialize lists and arrays to store data
# chiN = []
# r1ab = []
# chiN_data = {}
# with open(sys.argv[1], "r") as file:
#     # Step 1: Read the header line
#     header = file.readline().strip().split()
    
#     # The first column is 'r1ab', the rest are 'chiN' values
#     chiN = [int(col.replace('chiN', '')) for col in header[1:]]
    
#     # Step 2: Read the data lines
#     data = np.loadtxt(file)
    
#     # Separate the data into r1ab and chiN arrays
#     r = data[:, 0]
#     for i, chiN_value in enumerate(chiN):
#         chiN_data[f'chiN{chiN_value}'] = data[:, i+1]

# # Step 3: Assign chiN arrays to individual variables
# hr1 = chiN_data[f'chiN{chiN[0]}']
# hr2 = chiN_data[f'chiN{chiN[1]}']
# hr3 = chiN_data[f'chiN{chiN[2]}']


# print(hr3)
chiN = [17, 33, 48]

params1 = [0.8,	4.3, 4.8, 3.2]
params2 = [1.5, 4.3, 4.8, 3.2]
params3 = [2.5, 4.3, 4.8, 3.2]
# params2 = [0.7,	3.6, 3.4, 5.75210487]
# params3 = [1.15, 3.6, 3.4, 5.75210487]
# # # params5 = [2.45, 1.46, 4.1, 3.7]

cmap1 = plt.get_cmap('CMRmap')
# params_skew = [params1, params2, params3, params3]

# h_r=[hr1, hr2, hr3]
# Define Gaussian function
# def gaussian(r, A, r_peak, sigma):
#     return A * np.exp(-((r - r_peak)**2) / (2 * sigma**2))

# # Define Lorentzian function
# def lorentzian(r, A, r_peak, gamma):
#     return A / (1 + ((r - r_peak) / gamma)**2)

# Define Skew-normal function
def skew_normal(r, A, r_peak, sigma, alpha):
    return A * skewnorm.pdf(r, alpha, loc=r_peak, scale=sigma)

# # Define Generalized Normal function
# def generalized_normal(r, A, r_peak, beta, p):
#     return A * np.exp(-np.abs((r - r_peak) / beta)**p)

# for i in range(3):
# Fit Gaussian
# params_gauss, _ = curve_fit(gaussian, r, h_r[0])
# A_gauss, r_peak_gauss, sigma_gauss = params_gauss

    # # Fit Lorentzian
    # params_lorentz, _ = curve_fit(lorentzian, r, h_r[i])
    # A_lorentz, r_peak_lorentz, gamma_lorentz = params_lorentz

    #Fit Skew-normal
    # params_skew, _ = curve_fit(skew_normal, r, h_r[i], p0=[1, 5, 1, 0])
    # A_skew, r_peak_skew, sigma_skew, alpha_skew = params_skew
    # print(params_skew)

    # # Fit Generalized Normal
    # params_gen_norm, _ = curve_fit(generalized_normal, r, h_r[i], p0=[1, 5, 1, 2])
    # A_gen_norm, r_peak_gen_norm, beta_gen_norm, p_gen_norm = params_gen_norm

    # Plot data and fits
# params_skew = [2.2, 1.94744477, 4.31078809, 4.8]
# plt.plot(r, hr1, 'o', label=r'$\chi$N=%d'%chiN[i], alpha=0.5, markersize=2, color=cmap1(0.23+float(i)/4))
plt.plot(r, hr1, alpha=0.5,color=cmap1(0.23+float(0)/4))
plt.plot(r, hr2, alpha=0.5, color=cmap1(0.23+float(1)/4))
plt.plot(r, hr3, alpha=0.5,color=cmap1(0.23+float(2)/4))
# plt.plot(r, gaussian(r, *params_gauss), '-', label='Gaussian fit')
# plt.plot(r, lorentzian(r, *params_lorentz), '--', label='Lorentzian fit')
plt.plot(r, skew_normal(r, *params1), '-.', color=cmap1(0.23+float(0)/4))
plt.plot(r, skew_normal(r, *params2), '-.', color=cmap1(0.23+float(1)/4))
plt.plot(r, skew_normal(r, *params3), '-.', color=cmap1(0.23+float(2)/4))
plt.plot()
    # plt.plot(r, generalized_normal(r, *params_gen_norm), ':', label='Generalized Normal fit')
# np.savetxt("skew_fit.txt", np.column_stack((r, skew_normal(r, *params_skew[0]), skew_normal(r,*params_skew[1]), skew_normal(r, *params_skew[2]))), fmt='%.5f')
# plt.xlim(0, 20)
plt.ylim(0,0.4)
plt.xlabel('r')
plt.ylabel('h(r)')
plt.legend()
plt.savefig('skew_fit.pdf')
plt.show()

header1 = f"r1ab chiN{chiN[0]} chiN{chiN[1]}"
header2 = f"{params1}\n{params2}\n"
with open("skew_fit.txt", "w") as f:
    f.write(header1)
    f.write("\n")
    f.write(header2)
    f.write("\n")
    np.savetxt(f, np.column_stack((r, skew_normal(r, *params1), skew_normal(r,*params2), skew_normal(r,*params3))), fmt='%.5f')
    f.close()

