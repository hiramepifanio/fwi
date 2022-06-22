import numpy as np
import matplotlib.pyplot as plt
import sys

script = sys.argv[1]
folderPath = "../data/" + script + "/"

file_name = sys.argv[2]
data = np.fromfile(folderPath + file_name+".bin", 'float32')
shape = np.loadtxt(folderPath + file_name+".config", dtype=int)
data = data.reshape(shape)

#plt.plot(data)
plt.imshow(data)
#plt.imshow(data, aspect='auto')
#plt.clim(-0.0001,0.0001)
#plt.clim(-data.max()/50, data.max()/50)
plt.clim(-10,10)
plt.colorbar()
plt.show()