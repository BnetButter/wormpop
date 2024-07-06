import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('summary.tsv', delimiter='\t')

X = np.array(data["Timestep"])
Y = np.array(data["Number Worms"])

# Baseline correction
# Assuming the last 10% of the dataset is where the signal should have decayed
baseline_value = np.mean(Y[int(0.9 * len(Y)):])
corrected_fid = Y - baseline_value

Y_fft = np.fft.fft(Y)
# Compute the frequencies
n = len(corrected_fid)          # Length of the signal
sample_rate = 1     # Sampling frequency in Hz
freq = np.fft.fftfreq(n, d=1)

# Plot the magnitudes of the FFT
plt.figure()
plt.stem(freq, np.abs(Y_fft), 'b', markerfmt=" ", basefmt="-b")
plt.title('FFT of the signal')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude')
plt.grid()
plt.legend()
plt.grid(True)
plt.savefig('windowed_fid.png')