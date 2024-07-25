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
# Constants


def freqMix(real_data, num_samples):
    samples_to_return = []
    for i in range(num_samples):
        # iterate over each frequency and pick one
        new_signal = []
        last_Selection = None
        for frequency_sample in range(len(real_data[0])):

            if last_Selection != None:
                filteredList = filter_list(
                    real_data[:, frequency_sample], last_Selection-10, last_Selection+10)
                try:
                    selection = random.choice(filteredList)
                except:
                    selection = random.choice(real_data[:, frequency_sample])

            else:
                selection = random.choice(real_data[:, frequency_sample])
            new_signal.append(selection)
            last_Selection = selection
        samples_to_return.append(new_signal)
    return samples_to_return


def filter_list(list, lower_bound, upper_bound):
    return_list = []
    for item in list:
        if item < upper_bound and item > lower_bound:
            return_list.append(item)
    return return_list


power_options = [238]  # [53, 62, 74, 87, 103, 107,
# 110, 114, 122, 135, 144, 170, 201, 238]
for lasing_power in power_options:
    # lasing_power = 62
    PULSE_FILENAME = "C:/Users/Sean/Documents/Programing/Python/PURDUE/RealData/RawData" + \
        str(lasing_power)+"[uJcm2].csv"
    FREQUENCIES_FILENAME = "C:/Users/Sean/Documents/Programing/Python/PURDUE/RealData/Frequencies.csv"
    dataset = read_csv(PULSE_FILENAME, header=None)
    frequencies = read_csv(FREQUENCIES_FILENAME, header=None)
    # Convert dataset into the format rows:pulse instead of collumn:pulse
    time_signals = dataset.T.reset_index().reindex().values[:, 1:]
    # print(time_signals)
    frequencies = frequencies.values[0, :]
    fig, axs = plt.subplots(4)
    # plt.margins(0, None)
    fig.subplots_adjust(top=0.94)
    # fig.tight_layout()  # h_pad=23)
    # Graph an example Power Spectrum
    plt.figure()
    # plt.plot(frequencies, time_signals[4])
    # plt.plot(frequencies, time_signals[2])
    # plt.plot(frequencies, time_signals[0])
    # plt.title("Real Data")
    axs[0].matshow(time_signals)
    axs[0].set_xlabel('Frequency (Hz)')
    axs[0].set_ylabel('Sample')
    axs[0].set_yticks(np.arange(0, 100.1, 100/2))
    axs[0].set_title('Real Data: ' + str(lasing_power) + "[uJcm2]")
    axs[0].grid(True)

    # # Define Direchlet ethod

    def DirechletLIM(num_to_gen, fake_data, num_to_pull):
        size = len(fake_data)
        signals_to_return = []

        signals = [random.randrange(0, size) for i in range(num_to_pull)]
        b = [0.1*((i+1)**(1/5))
             for i in range(len(signals))]  # 0.1 #0.1*((i+1)**(1/5))
        alpha = [random.choice(b) for i in range(len(signals))]
        for j in range(num_to_gen):
            distributions = dirichlet.rvs(alpha, size=1)
            new_signal = [0 for i in range(len(frequencies))]
            c = 0
            for i in signals:
                new_signal += np.array(fake_data[i]) * distributions[0][c]
                c += 1
            signals_to_return.append(new_signal)
        return signals_to_return
    # Define average + Normal

    def Direchlet(num_to_gen, fake_data):
        signals_to_return = []
        alpha = [1 for i in range(len(fake_data))]
        for j in range(num_to_gen):
            distributions = dirichlet.rvs(alpha, size=1)
            new_signal = [0 for i in range(len(frequencies))]
            for i in range(0, len(fake_data)):
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
        for j in range(num_to_gen):
            new_signal = []
            for i in range(len(average)):
                new_signal.append(
                    (average[i] + random.gauss(0.0, average_deviation)))
            signals_to_return.append(new_signal)
            # print(len(new_signal))
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

    def Dists(num_to_gen, fake_data):
        dist_list = []
        signals_to_return = []
        for wavelength in range(len(fake_data[0])):
            WLnums = []
            for pulse in range(len(fake_data)):
                WLnums.append(fake_data[pulse][wavelength])
            dist = stats.norm  # fold norm
            # [100,200][0,0.75] kinda worked with [100,400]
            dist = stats.fit(dist, WLnums, [(0, 400), (0, 0.75)])
            # dist.plot()
            # plt.show()
            dist_list.append(dist)
        for i in range(num_to_gen):
            signal = []
            for dist in dist_list:
                wl = stats.norm.rvs(
                    loc=dist.params.loc, scale=dist.params.scale, size=1, random_state=None)
                signal.append(wl)
            signals_to_return.append(signal)
        return signals_to_return

    def mse(ground_truth, gen_data):
        # loop through each pair of shots
        meanmse = 0
        for shot in range(len(gen_data)):
            mean = 0
            for amp in range(len(gen_data[shot])):
                error = ground_truth[shot][amp] - gen_data[shot][amp]
                serror = error ** 2
                mean += serror
            mean /= len(gen_data[shot])
            meanmse += mean
        return meanmse/len(gen_data)

        # calculate MSE between them
        # average all MSE's

    # np.savetxt("110done.csv", time_signals, delimiter=",")
    gen_samples = DirechletLIM(100, time_signals, 5)
    gen_samples = np.rot90(np.array(gen_samples), 3)
    np.savetxt('direchlet.csv', gen_samples, delimiter=",")
    # print("Direchlet MSE:" + str(lasing_power) +
    #       '[ujcm2]'+str(mse(time_signals, gen_samples)))
    # gen_samples = Direchlet(100, time_signals)
    # plt.plot(frequencies, gen_samples[0])
    # plt.plot(frequencies, gen_samples[1])
    # plt.plot(frequencies, gen_samples[2])
    # plt.plot(frequencies, time_signals[4], color='b', label="Real Data")
    # # plt.plot(frequencies, time_signals[2], color='b')
    # # plt.plot(frequencies, time_signals[0], color='b')
    # plt.plot(frequencies, gen_samples[0], color='r', label="Direchlet Data")
    # # plt.plot(frequencies, gen_samples[1], color='r')
    # # plt.plot(frequencies, gen_samples[2], color='r')
    # plt.xlabel("Wavelength")
    # plt.ylabel("Power")
    # plt.title("Overlay Generated and Real Data " +
    #           str(lasing_power)+"[uJcm2]")
    # plt.legend(loc="upper left")
    axs[1].matshow(gen_samples)
    axs[1].set_xlabel('Frequency Sample')
    axs[1].set_ylabel('Shots')
    axs[1].set_title('Direchlet Data')
    axs[1].set_yticks(np.arange(0, 100.1, 100/2))
    axs[1].grid(True)

    gen_samples = averageWithDist(100, time_signals, 10000, False)
    gen_samples = np.rot90(np.array(gen_samples), 3)
    np.savetxt('avrgwDist.csv', gen_samples, delimiter=",")

    plt.plot(frequencies, time_signals[4], color='b', label="Real Data")
    # plt.plot(frequencies, time_signals[2], color='b')
    # plt.plot(frequencies, time_signals[0], color='b')
    plt.plot(frequencies, gen_samples[0], color='r', label="Direchlet Data")
    # plt.plot(frequencies, gen_samples[1], color='r')
    # plt.plot(frequencies, gen_samples[2], color='r')
    plt.xlabel("Wavelength")
    plt.ylabel("Power")
    plt.title("Overlay Generated and Real Data " +
              str(lasing_power)+"[uJcm2]")
    plt.legend(loc="upper left")

    print("Averge WIth DIST ist MSE:"+str(lasing_power) +
          '[ujcm2]'+str(mse(time_signals, gen_samples)))
    axs[2].matshow(gen_samples)
    axs[2].set_xlabel('Frequency Sample')
    axs[2].set_ylabel('Shots')
    axs[2].set_yticks(np.arange(0, 100.1, 100/2))

    axs[2].set_title('Average Dist Data')
    axs[2].grid(True)

    # gen_samples = averageWithNormal(100, time_signals)
    # # plt.plot(frequencies, gen_samples[0])
    # # plt.plot(frequencies, gen_samples[1])
    # # plt.plot(frequencies, gen_samples[2])
    # axs[3].matshow(gen_samples)
    # axs[3].set_xlabel('Frequency Sample')
    # axs[3].set_ylabel('Shots')
    # axs[3].set_yticks(np.arange(0, 100.1, 100/2))
    # axs[3].set_title('Average Norm Data')
    # axs[3].grid(True)

    gen_samples = freqMix(time_signals, 100)
    print("Freq Mix MSE:"+str(lasing_power) +
          '[ujcm2]'+str(mse(time_signals, gen_samples)))
    # gen_samples = Dists(100, time_signals)
    # plt.plot(frequencies, gen_samples[0])
    # plt.plot(frequencies, gen_samples[1])
    # plt.plot(frequencies, gen_samples[2])
    axs[3].matshow(gen_samples)
    axs[3].set_xlabel('Frequency Sample')
    axs[3].set_ylabel('Shots')
    axs[3].set_yticks(np.arange(0, 100.1, 100/2))
    axs[3].set_title('Wavelength Mix Data')
    axs[3].grid(True)
    plt.show()
    # fig.savefig('C:/Users/Sean/Documents/Programing/Python/PURDUE/IMAGES/Complist' +
    # str(lasing_p150ower) + '[ujcm2].png')
