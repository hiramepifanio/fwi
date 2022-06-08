import numpy as np
import matplotlib.pyplot as plt
import sys

script = sys.argv[1]
folderPath = "../data/" + script + "/"

data = np.fromfile(folderPath + "grad.bin", 'float32')
shape = np.loadtxt(folderPath + "grad.config", dtype=int)
data = data.reshape(shape[::-1])
data = data.T

plt.imshow(data)
plt.clim(-0.0001,0.0001)
#plt.clim(-2,2)
plt.colorbar()
plt.show()