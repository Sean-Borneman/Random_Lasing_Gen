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
time_signals = []
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
    time_signals.append(unique_signal)
# Graph an example Power Spectrum
plt.figure(figsize=(10, 5))
plt.plot(frequencies, time_signals[4])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power')
plt.title('Example Fake Power Spectrum')
plt.grid(True)
plt.xlim(0)
# Define Direchlet ethod


def Direchlet(num_to_gen, fake_data, INITIAL_NUM_FAKE_SAMPLES):
    signals_to_return = []
    alpha = [1 for i in range(INITIAL_NUM_FAKE_SAMPLES)]
    for j in range(num_to_gen):
        distributions = dirichlet.rvs(alpha, size=1)
        new_signal = [0 for i in range(int(1340))]#total_time*sampling_rate
        for i in range(0, INITIAL_NUM_FAKE_SAMPLES):
            new_signal += np.array(fake_data[i]) * distributions[0][i]
        signals_to_return.append(new_signal)
    return signals_to_return
# Define average + Normal


def averageWithNormal(num_to_gen, fake_data):
    signals_to_return = []
    total_deviation = 0
    sum = np.sum(fake_data, axis=0)
    average = sum/len(fake_data)

    naverage = np.array(average)
    ntime_sig = np.array(average)
    for i in range(len(fake_data)):
        deviation = naverage - ntime_sig[i]
        # Returns differant output than original function VERY SUS
        # TODO:PLEASE FIX DIFFERANT OUTPUTS
        total_deviation += np.sum(np.absolute(deviation))

    average_deviation = total_deviation / (len(ntime_sig)*len(ntime_sig))
    new_signal = []
    for j in range(num_to_gen):
        for i in range(len(average)):
            new_signal.append(
                average[i] + random.gauss(0.0, average_deviation))
        signals_to_return.append(new_signal)
    return signals_to_return
# Define Average + fitdist


def averageWithDist(num_to_gen, fake_data, DEVIATION_THRESHOLD, DEVIATION_FILTER):
    signals_to_return = []
    deviations = []
    counter = 0
    total_deviation = 0
    sum = 0
    for i in range(len(fake_data)):
        sum += fake_data[i]
    average = sum/len(fake_data)
    for i in range(len(fake_data)):
        for j in range(len(fake_data[i])):
            counter += 1
            deviation = average[j] - fake_data[i][j]
            total_deviation += math.sqrt((average[j] - fake_data[i][j])**2)
            if ((abs(deviation) > DEVIATION_THRESHOLD) and (DEVIATION_FILTER == True)):
                pass
            else:
                deviations.append(deviation)
    dist = stats.norm
    bounds = [(-30, 30), (0, 100)]
    # bounds = [(0,30), (0,1)]
    dist = stats.fit(dist, deviations, bounds)

    for j in range(num_to_gen):
        new_signal = []
        for i in range(len(average)):
            sample = stats.norm.rvs(
                loc=dist.params.loc, scale=dist.params.scale, size=1, random_state=None)
            # sample = stats.dgamma.rvs(dist.params.a, dist.params.loc, dist.params.scale, size=1, random_state=None)
            new_signal.append(average[i] + sample[0])
        signals_to_return.append(new_signal)
    return signals_to_return


# Time Each method
# size of num_to_gen to see scaling
test_sizes = [10, 100, 1000, 10000, 100000]

dirichlet_times = []
averageNorm_times = []
averageDist_times = []
gen_sample = Direchlet(100, time_signals, INITIAL_NUM_FAKE_SAMPLES)
gen_sample = np.rot90(np.array(gen_sample), 3)
np.savetxt('Direchlet.csv', gen_sample, delimiter=",")

gen_sample = averageWithNormal(100, time_signals)
gen_sample = np.rot90(np.array(gen_sample), 3)
np.savetxt('averageNorm.csv', gen_sample, delimiter=",")
# for size in test_sizes:
#     start_time = t.time()
#     Direchlet(size, time_signals, INITIAL_NUM_FAKE_SAMPLES)
#     end_time = t.time()
#     dirichlet_times.append(end_time-start_time)

#     start_time = t.time()
#     averageWithNormal(size, time_signals)
#     end_time = t.time()
#     averageNorm_times.append(end_time-start_time)
#     # start_time = t.time()
#     # averageWithDist(size, time_signals, DEVIATION_THRESHOLD, DEVIATION_FILTER)
#     # end_time = t.time()
#     # averageNorm_times.append(end_time-start_time)
# print(dirichlet_times)
# print(averageNorm_times)
