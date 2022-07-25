import numpy as np
import matplotlib.pyplot as plt

file_name = 'gradr'
folderPath = "../data/"

grad = np.loadtxt(folderPath + file_name+'.data', dtype=np.float32)
grad_shape = np.loadtxt(folderPath + file_name+'.shp', dtype=int)
grad = grad.reshape(grad_shape)

figname = 'grads-'+file_name
fig = plt.figure(figname)
plt.imshow(grad, aspect='auto')

plt.xlabel('$j$')
plt.ylabel('$i$')
plt.clim(-7e9,7e9)
cbar = plt.colorbar(mappable = None)
cbar.set_label('$(\\nabla\\Phi)_{ij}\\ \\frac{N^2s}{m^5}$')

plt.show()
fig.savefig(figname)