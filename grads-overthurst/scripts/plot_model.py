import numpy as np
import matplotlib.pyplot as plt

file_name = 'vel'
folderPath = "../data/"

vel = np.loadtxt(folderPath + file_name+'.data', dtype=np.float32)
vel_shape = np.loadtxt(folderPath + file_name+'.shp', dtype=int)
vel = vel.reshape(vel_shape)

vel0 = np.loadtxt(folderPath + file_name+'0.data', dtype=np.float32)
vel0 = vel0.reshape(vel_shape)

sou = np.loadtxt(folderPath + 'source.data', dtype=np.float32)
shape = np.loadtxt(folderPath + 'source.shp', dtype=int)
sou = sou.reshape(shape)

file_name = 'recs'
recs = np.loadtxt(folderPath + file_name+'.data', dtype=np.float32)
recs_shape = np.loadtxt(folderPath + file_name+'.shp', dtype=int)
recs = recs.reshape(recs_shape)

figname = 'grads-vel'
fig = plt.figure(figname, figsize=(15, 5))

plt.subplot(121)
plt.title('Mapa de velocidade real')
plt.imshow(vel, aspect='auto')
plt.plot(recs[1], recs[0],'vg', label = 'receptores')
plt.plot(sou[1], sou[0],'*r', label = 'fonte')
plt.xlabel('$x \\ (\\times 25\\ m$)')
plt.ylabel('$z \\ (\\times 25\\ m$)')
plt.legend()
plt.clim(3000, 6000)
cbar = plt.colorbar(mappable = None)
cbar.set_label('$v\\ (m/s)$')


plt.subplot(122)
plt.title('Modelo inicial')
plt.imshow(vel0, aspect='auto')
plt.plot(recs[1], recs[0],'vg', label = 'receptores')
plt.plot(sou[1], sou[0],'*r', label = 'fonte')
plt.xlabel('$x \\ (\\times 25\\ m$)')
plt.ylabel('$z \\ (\\times 25\\ m$)')
plt.legend()
plt.clim(3000, 6000)
cbar = plt.colorbar(mappable = None)
cbar.set_label('$v\\ (m/s)$')


fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=None)

plt.show()
fig.savefig(figname)