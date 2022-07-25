import numpy as np
import matplotlib.pyplot as plt

file_name = 'wave_'
folder_path = "../data/"

sou = np.loadtxt(folder_path + 'source.data', dtype=np.float32)
shape = np.loadtxt(folder_path + 'source.shp', dtype=int)
sou = sou.reshape(shape)

figname = 'sp-waves'
step = 60
fig = plt.figure(figname, figsize=(12, 15))
for i in range(0, 6):
    wave = np.loadtxt(folder_path + file_name + str(step*i) + '.data', dtype=np.float32)
    shape = np.loadtxt(folder_path + file_name + str(step*i) + '.shp', dtype=int)
    wave = wave.reshape(shape)

    plt.subplot(321 + i)
    plt.imshow(wave, aspect='auto')
    plt.plot(sou[1], sou[0],'*r')
    plt.xlabel('$x \\ (\\times 15\\ m$)')
    plt.ylabel('$z \\ (\\times 15\\ m$)')
    plt.clim(-10,10)
    legend = 't = ' + str(step*i*3) + ' ms'
    plt.text(6, 15, legend, bbox={'facecolor': 'white', 'pad': 5})
    cbar = plt.colorbar(mappable = None)
    cbar.set_label('$p\\ (N/m^2)$')

fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=0.4)

#plt.show()
fig.savefig(figname)