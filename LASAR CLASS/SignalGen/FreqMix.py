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
from pandas import read_csv

lasing_power = 238
PULSE_FILENAME = "C:/Users/Sean/Documents/Programing/Python/PURDUE/RealData/RawData" + \
    str(lasing_power)+"[uJcm2].csv"
FREQUENCIES_FILENAME = "C:/Users/Sean/Documents/Programing/Python/PURDUE/RealData/Frequencies.csv"
dataset = read_csv(PULSE_FILENAME, header=None)
frequencies = read_csv(FREQUENCIES_FILENAME, header=None)
# Convert dataset into the format rows:pulse instead of collumn:pulse
shots = dataset.T.reset_index().reindex().values[:, 1:]

frequencies = frequencies.values[0, :]

# fig, axs = plt.subplots(5)
# Graph an example Power Spectrum
# plt.figure(figsize=(10, 5))
# plt.plot(frequencies, time_signals[4])
# plt.plot(frequencies, time_signals[2])
# plt.plot(frequencies, time_signals[0])
# plt.title("Real Data")
# axs[0].matshow(time_signals)
# axs[0].set_xlabel('Frequency (Hz)')
# axs[0].set_ylabel('Sample')
# axs[0].set_yticks(np.arange(0, 100.1, 100/2))
# axs[0].set_title('Real Data: ' + str(lasing_power) + "[uJcm2]")
# axs[0].grid(True)


def freqMix(real_data, num_samples):
    new_dataset = []
    for i in range(num_samples):
        # iterate over each frequency and pick one
        new_signal = []
        last_Selection = None
        for frequency_sample in range(len(real_data[0])):

            if last_Selection != None:
                filteredList = filter_list(
                    real_data[:, frequency_sample], last_Selection-1000, last_Selection+1000)
                selection = random.choice(filteredList)
            else:
                selection = random.choice(real_data[:, frequency_sample])
            new_signal.append(selection)
            last_Selection = selection
        new_dataset.append(new_signal)
    return new_dataset


def filter_list(list, lower_bound, upper_bound):
    return_list = []
    for item in list:
        if item < upper_bound and item > lower_bound:
            return_list.append(item)
    return return_list


gen_sample = freqMix(shots, 100)


np.savetxt('freqMIX.csv', gen_sample, delimiter=",")
plt.figure(figsize=(10, 5))
plt.plot(shots[0], color='b', label="Real Data")
plt.plot(gen_sample[0], color='r', label="Generated Data")

plt.xlabel('Frequency Num')
plt.ylabel('Aplitude')
plt.title('Mix pulse')
plt.grid(True)
plt.legend(loc="upper left")
plt.xlim(0)
plt.show()
gen_sample = np.rot90(np.array(gen_sample), 3)
np.savetxt('freqMIX.csv', gen_sample, delimiter=",")
# import numpy as np
# ground_truth = [[0, 1, 1, 0], [2, 3, 4, 5]]
# gen_data = [[0, 0, 2, 3], [9, 6, 7, 1]]
# a = np.array(ground_truth)  # your x
# b = np.array(gen_data)  # your y9
# print((a-b)**2).mean(axis=1)
