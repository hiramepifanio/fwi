import numpy as np
import matplotlib.pyplot as plt

file_name = 'velb'
folderPath = "../data/"

velb = np.loadtxt(folderPath + file_name+'.data', dtype=np.float32)
velb_shape = np.loadtxt(folderPath + file_name+'.shp', dtype=int)
velb = velb.reshape(velb_shape)

sou = np.loadtxt(folderPath + 'source.data', dtype=np.float32)
shape = np.loadtxt(folderPath + 'source.shp', dtype=int)
sou = sou.reshape(shape)

borda = 50
nz = velb_shape[0] - 2*borda
nx = velb_shape[1] - 2*borda
ret_x=[borda,borda+nx,borda+nx,borda,borda]
ret_z=[borda,borda,borda+nz,borda+nz,borda]

figname = 'inf-velb'
fig = plt.figure(figname)
plt.imshow(velb, aspect='auto')
plt.plot(sou[1] + borda, sou[0] + borda,'*r', label = 'fonte')
plt.plot(ret_x,ret_z,'r', label = 'borda')
plt.xlabel('$x \\ (\\times 15\\ m$)')
plt.ylabel('$z \\ (\\times 15\\ m$)')
plt.legend()
plt.clim(2000, 3000)
cbar = plt.colorbar(mappable = None)
cbar.set_label('$v\\ (m/s)$')
#plt.show()
fig.savefig(figname)