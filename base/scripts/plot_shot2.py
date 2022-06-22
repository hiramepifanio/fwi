import numpy as np
import matplotlib.pyplot as plt
import sys

script = sys.argv[1]
folderPath = "../data/" + script + "/"

datad = np.fromfile(folderPath + "shotd"+".bin", 'float32')
datar = np.fromfile(folderPath + "shotr"+".bin", 'float32')
shape = np.loadtxt(folderPath + "shotd"+".config", dtype=int)

datad = datad.reshape(shape[::-1])
datad = datad.T

datar = datar.reshape(shape[::-1])
datar = datar.T

print("shots")
print((datad**2).mean())
print((datar**2).mean())
print(((datad-datar)**2).mean())

lim = 0.03
fig = plt.figure()
plt.imshow(datad,aspect='auto',cmap='binary');plt.clim(-lim, lim)
plt.title('Campo direto')
fig.savefig('shotd')

fig = plt.figure()
plt.imshow(datar,aspect='auto',cmap='binary');plt.clim(-lim, lim)
plt.title('Campo reconstruido')
fig.savefig('shotr')

fig = plt.figure()
plt.imshow(datad-datar,aspect='auto',cmap='binary');plt.clim(-lim, lim)
plt.title('Diferença')
fig.savefig('shotdiff')