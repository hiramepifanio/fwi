import numpy as np
import matplotlib.pyplot as plt
import sys

script = sys.argv[1]
folderPath = "../data/" + script + "/"

datad = np.fromfile(folderPath + "grad"+".bin", 'float32')
datar = np.fromfile(folderPath + "gradr"+".bin", 'float32')
shape = np.loadtxt(folderPath + "grad"+".config", dtype=int)

datad = datad.reshape(shape)

datar = datar.reshape(shape)

lim = 100
fig = plt.figure()
plt.imshow(datad,aspect='auto',cmap='binary');plt.clim(-lim, lim)
plt.title('Campo direto')
fig.savefig('grad')

fig = plt.figure()
plt.imshow(datar,aspect='auto',cmap='binary');plt.clim(-lim, lim)
plt.title('Campo reconstruido')
fig.savefig('gradr')

fig = plt.figure()
plt.imshow(datad-datar,aspect='auto',cmap='binary');plt.clim(-lim, lim)
plt.title('Diferença')
fig.savefig('graddiff')