import numpy as np
import matplotlib.pyplot as plt
import sys

script = sys.argv[1]
folderPath = "../data/" + script + "/"

file_name = "grad"
data = np.fromfile(folderPath + file_name+".bin", 'float32')
shape = np.loadtxt(folderPath + file_name+".config", dtype=int)
data = data.reshape(shape)

#plt.plot(data)
plt.imshow(data)
plt.clim(-0.0001,0.0001)
#plt.clim(-2,2)
plt.colorbar()
plt.show()