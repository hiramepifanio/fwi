import numpy as np
import matplotlib.pyplot as plt
import sys

file_name = 'ricker_wavelet'
folderPath = "../data/"

data = np.loadtxt(folderPath + file_name+".txt", dtype=np.float32)
shape = np.loadtxt(folderPath + file_name+".config", dtype=int)
print(shape)
data = data.reshape(shape)


x = np.arange(0, shape, 1)
fp = 13
t0 = 6/(np.pi*fp*2**0.5)
print(t0)

fig = plt.figure('ricker-wavelet')
plt.plot(x*2, data, label = '$f_p = 13$ Hz')
plt.xlabel('$t$ (ms)')
plt.ylabel('$s_0(t)$')
plt.legend()
plt.grid()
plt.show()
fig.savefig('ricker-wavelet')


dt = 0.002
freq=np.fft.rfftfreq(data.size,dt)
A = np.abs(np.fft.rfft(data))

fig = plt.figure('ricker-wavelet-fft')
plt.plot(freq, A)
plt.plot([fp, fp], [0,16], '--', label = '$f_p$')
plt.plot([3*fp, 3*fp], [0,16], '--', label = '$3f_p$')
plt.xlabel('$f$ (Hz)')
plt.ylabel('$c_f$')
plt.xlim(0,50)
plt.legend()
plt.grid()
plt.show()
fig.savefig('ricker-wavelet-fft')