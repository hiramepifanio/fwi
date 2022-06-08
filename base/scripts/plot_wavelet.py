import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../data/sincint_wavelet.txt", "float")



print(data.shape)
plt.plot(data)
plt.show()
