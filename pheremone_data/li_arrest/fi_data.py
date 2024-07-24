import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# Load the data
filename = "arrest-data.csv"
data = pd.read_csv(filename)

# Extract the data
x = data['x']
y = data['Curve1']


def starve_from_l1_arrest(num_days):
    L = 98.09173354410315
    x0 = 17.484181571141384
    k = -0.47620124366404964

        # Define the reverse sigmoid function
    def reverse_sigmoid(x, L, x0, k):
        return L / (1 + np.exp(-k * (x - x0)))

    return reverse_sigmoid(num_days, L, x0, k)


# Initial guess for the parameters
initial_guess = [max(y), np.median(x), 1]

# Fit the curve
params, covariance = curve_fit(reverse_sigmoid, x, y, p0=initial_guess)

# Extract the parameters
L, x0, k = params

print(f"L: {L}")
print(f"x0: {x0}")
print(f"k: {k}")


# Generate the fitted curve
x_fit = np.linspace(min(x), max(x), 100)
y_fit = reverse_sigmoid(x_fit, L, x0, k)

# Plot the data and the fitted curve
plt.scatter(x, y, label='Data')
plt.plot(x_fit, y_fit, color='red', label='Fitted curve')
plt.xlabel('x')
plt.ylabel('Curve1')
plt.legend()
plt.title('Reverse Sigmoid Curve Fitting')
plt.savefig('reverse_sigmoid_fit.png')