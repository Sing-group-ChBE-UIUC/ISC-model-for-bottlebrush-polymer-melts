import numpy as np
import matplotlib.pyplot as plt
import sys
from sklearn.linear_model import RidgeCV

# Generate noisy data
np.random.seed(0)
x_noisy, y_noisy = np.loadtxt(sys.argv[1], unpack=True)
x_sample, y_sample = np.loadtxt(sys.argv[2], unpack=True)

# Fit a ridge regression model to smooth the noisy data while constraining it to pass through the sample data
model = RidgeCV(alphas=np.logspace(-100, 100, 10000))
model.fit(x_sample[:, np.newaxis], y_sample)

# Predict smoothed data
y_smoothed = model.predict(x_noisy[:, np.newaxis])

# Plot original noisy data, sample data, and smoothed data
plt.figure(figsize=(10, 6))
plt.plot(x_noisy, y_noisy, label='Noisy Data', color='gray')
plt.plot(x_sample, y_sample, 'o', label='Sample Data', color='red')
plt.plot(x_noisy, y_smoothed, label='Smoothed Data', color='blue')
plt.legend()
plt.title('Constrained Smoothing with Ridge Regression')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.show()
