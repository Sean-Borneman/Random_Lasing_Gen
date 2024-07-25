# Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import random
from scipy.stats import dirichlet
import math
import numpy as np
from scipy import stats
import time as t
# Constants
INITIAL_NUM_FAKE_SAMPLES = 100
FREQUENCY_ERROR = 0
AMPLITUDE_ERROR = 0
GENRAL_ERROR_AMPLITUDE = 10
DEVIATION_THRESHOLD = 100
DEVIATION_FILTER = True
total_time = 2.0      # Total time duration of the signal in seconds
sampling_rate = 1000
# Generate a set of Fake data (replace with real data once assecible)
time = np.linspace(0, total_time, int(
    sampling_rate * total_time), endpoint=False)
frequencies = np.fft.fftfreq(len(time), 1 / sampling_rate)
frequency_signals = []
for i in range(INITIAL_NUM_FAKE_SAMPLES):
    sampling_rate = 1000  # Sampling rate in Hz (samples per second)
    # Frequency of the signal in Hz
    classic_frequency = 5.0 + random.uniform(-FREQUENCY_ERROR, FREQUENCY_ERROR)
    # Amplitude of the signal
    classic_amplitude = 3.0 + random.uniform(-AMPLITUDE_ERROR, AMPLITUDE_ERROR)
    lasing_frequency = 10.0 + random.uniform(-FREQUENCY_ERROR, FREQUENCY_ERROR)
    lasing_amplitude = 5.0 + random.uniform(-AMPLITUDE_ERROR, AMPLITUDE_ERROR)

    # Generate a time-domain signal (sinusoidal signal)
    time_domain_signal = classic_amplitude * \
        np.sin(2 * np.pi * classic_frequency * time)
    time_domain_signal += lasing_amplitude * \
        np.sin(2 * np.pi * lasing_frequency * time)
    unique_signal = np.abs(np.fft.fft(time_domain_signal))
    for i in range(len(unique_signal)):
        # add uniform error
        unique_signal[i] = unique_signal[i] + \
            random.uniform(-1*GENRAL_ERROR_AMPLITUDE, GENRAL_ERROR_AMPLITUDE)

    if (unique_signal[0] > 10000):
        print("ALERT")
        print(classic_frequency)
        print(classic_amplitude)
        print(lasing_frequency)
        print(lasing_amplitude)
    frequency_signals.append(unique_signal)

time_signals = []
for signal in frequency_signals:
    time_signals.append(signal)
    # time_signals.append(np.fft.ifft(signal))
full_dataset = []

for i in range(len(time_signals[0])):
    dataset = []
    for signal in frequency_signals:
        dataset.append(signal[i])
    full_dataset.append(dataset)

# Create a figure instance
fig = plt.figure(1, figsize=(9, 6))

# Create an axes instance
ax = fig.add_subplot(111)

# Create the boxplot
print(len(full_dataset))
bp = ax.boxplot(full_dataset[0:100])

# Save the figure
fig.savefig('fig1.png', bbox_inches='tight')
