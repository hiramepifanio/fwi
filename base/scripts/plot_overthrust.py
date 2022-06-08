import numpy as np
import matplotlib.pyplot as plt

data = np.fromfile("../data/overthrust.bin", 'float32')

nz, nx = 187, 801
data = data.reshape((nx, nz))
data = data.T

plt.imshow(data)
plt.show()
