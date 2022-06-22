import numpy as np
import matplotlib.pyplot as plt
import sys

script = sys.argv[1]
folderPath = "../data/" + script + "/"

data = np.fromfile(folderPath + "shot.bin", 'float32')
shape = np.loadtxt(folderPath + "shot.config", dtype=int)
data = data.reshape(shape[::-1])
data = data.T

plt.imshow(data, aspect = 'auto', cmap = 'binary')
#plt.clim(-data.max()/50,data.max()/50)
plt.clim(-2,2)
plt.colorbar()
plt.show()
