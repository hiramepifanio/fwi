import numpy as np
import matplotlib.pyplot as plt

file_name = 'buffer'
folderPath = "../data/"

buffer = np.loadtxt(folderPath + file_name+'.data', dtype=np.float32)
buffer_shape = np.loadtxt(folderPath + file_name+'.shp', dtype=int)
buffer = buffer.reshape(buffer_shape)

borda = 50
nz = buffer_shape[0] - 2*borda
nx = buffer_shape[1] - 2*borda
ret_x=[borda,borda+nx,borda+nx,borda,borda]
ret_z=[borda,borda,borda+nz,borda+nz,borda]

figname = 'inf-buffer'
fig = plt.figure(figname)
plt.imshow(buffer, aspect='auto')
plt.plot(ret_x,ret_z,'r', label = 'borda')
plt.xlabel('$x \\ (\\times 15\\ m$)')
plt.ylabel('$z \\ (\\times 15\\ m$)')
plt.legend()
#plt.clim(2000, 3000)
cbar = plt.colorbar(mappable = None)
cbar.set_label('$A$')
#plt.show()
fig.savefig(figname)