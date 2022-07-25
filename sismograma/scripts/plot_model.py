import numpy as np
import matplotlib.pyplot as plt

file_name = 'vel'
folderPath = "../data/"

vel = np.loadtxt(folderPath + file_name+'.data', dtype=np.float32)
vel_shape = np.loadtxt(folderPath + file_name+'.shp', dtype=int)
vel = vel.reshape(vel_shape)

sou = np.loadtxt(folderPath + 'source.data', dtype=np.float32)
shape = np.loadtxt(folderPath + 'source.shp', dtype=int)
sou = sou.reshape(shape)

file_name = 'recs'
recs = np.loadtxt(folderPath + file_name+'.data', dtype=np.float32)
recs_shape = np.loadtxt(folderPath + file_name+'.shp', dtype=int)
recs = recs.reshape(recs_shape)

figname = 'inf-model'
fig = plt.figure(figname)
plt.imshow(vel, aspect='auto')
plt.plot(recs[1], recs[0],'vg', label = 'receptores')
plt.plot(sou[1], sou[0],'*r', label = 'fonte')
plt.xlabel('$x \\ (\\times 15\\ m$)')
plt.ylabel('$z \\ (\\times 15\\ m$)')
plt.legend()
plt.clim(2000, 3000)
cbar = plt.colorbar(mappable = None)
cbar.set_label('$v\\ (m/s)$')
#plt.show()
fig.savefig(figname)