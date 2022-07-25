import numpy as np
import matplotlib.pyplot as plt

file_name = 'velocity_map'
folderPath = "../data/"

vel = np.loadtxt(folderPath + 'velocity_map.data', dtype=np.float32)
shape = np.loadtxt(folderPath + 'velocity_map.shp', dtype=int)
vel = vel.reshape(shape)

sou = np.loadtxt(folderPath + 'source.data', dtype=np.float32)
shape = np.loadtxt(folderPath + 'source.shp', dtype=int)
sou = sou.reshape(shape)
print(vel.min())
print(vel.max())

figname = 'sp-velocity-map'
fig = plt.figure(figname)
plt.imshow(vel, aspect='auto')
plt.plot(sou[1], sou[0],'*r', label = 'fonte')
plt.xlabel('$x \\ (\\times 15\\ m$)')
plt.ylabel('$z \\ (\\times 15\\ m$)')
plt.legend()
plt.clim(2000, 3000)
cbar = plt.colorbar(mappable = None)
cbar.set_label('$v\\ (m/s)$')
#plt.show()
fig.savefig(figname)