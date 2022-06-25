import numpy as np
import matplotlib.pyplot as plt
import sys

script = sys.argv[1]
folderPath = "../data/" + script + "/"

file_name = sys.argv[2]
data = np.fromfile(folderPath + file_name+".bin", 'float32')
shape = np.loadtxt(folderPath + file_name+".config", dtype=int)
data = data.reshape(shape)

plt.imshow(data, aspect = 'auto', cmap = 'binary')
plt.colorbar()
plt.show()