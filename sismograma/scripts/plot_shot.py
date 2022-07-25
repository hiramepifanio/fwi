import numpy as np
import matplotlib.pyplot as plt

file_name = 'shot'
folderPath = "../data/"

shot = np.loadtxt(folderPath + file_name+'.data', dtype=np.float32)
shot_shape = np.loadtxt(folderPath + file_name+'.shp', dtype=int)
shot = shot.reshape(shot_shape)

figname = 'inf-shot'
fig = plt.figure(figname)
plt.imshow(shot, aspect='auto', extent=[0,shot.shape[1],shot.shape[0]*3,0])
plt.xlabel('receptores')
plt.ylabel('$t\\ (ms)$')
plt.clim(-10, 10)
cbar = plt.colorbar(mappable = None)
cbar.set_label('$p\\ (N/m^2)$')
plt.show()
fig.savefig(figname)