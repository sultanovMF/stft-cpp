import numpy as np
import matplotlib.pyplot as plt

from koval_stft import stft, istft,  get_time_samples,  get_frequency_samples

signal = np.loadtxt("./signal_example.txt")

fs = 500
nperseg = 256
overlap = 128

spectogram = stft(signal, nperseg, overlap)


freq_samp = get_frequency_samples(nperseg // 2 + 1, fs) #only positive freq
time_samp = get_time_samples(spectogram.shape[1], overlap, fs)

plt.pcolormesh(time_samp, freq_samp, np.abs(spectogram), shading='gouraud')
plt.show()

