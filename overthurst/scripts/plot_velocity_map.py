import numpy as np
import matplotlib.pyplot as plt

file_name = 'vel'
folderPath = "../data/"

vel = np.loadtxt(folderPath + file_name+'.data', dtype=np.float32)
shape = np.loadtxt(folderPath + file_name+'.shp', dtype=int)
vel = vel.reshape(shape)

figname = 'over-vel'
fig = plt.figure(figname, figsize=(15,5))
plt.imshow(vel)
plt.xlabel('$x \\ (\\times 25\\ m$)')
plt.ylabel('$z \\ (\\times 25\\ m$)')
cbar = plt.colorbar(mappable = None)
cbar.set_label('$v\\ (m/s)$')
plt.show()
fig.savefig(figname)