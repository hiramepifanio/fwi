import numpy as np
import matplotlib.pyplot as plt
import sys

script = sys.argv[1]
folderPath = "../data/" + script + "/"

datad = np.fromfile(folderPath + "fieldd"+".bin", 'float32')
datar = np.fromfile(folderPath + "fieldr"+".bin", 'float32')
shape = np.loadtxt(folderPath + "fieldd"+".config", dtype=int)
datad = datad.reshape(shape)
datar = datar.reshape(shape)

#plt.plot(data)
#plt.imshow(datad)
# plt.imshow(datad, aspect='auto')
# plt.clim(-0.03,0.03)
# plt.title("Direct field at t = 100")
# #plt.clim(-data.max()/50, data.max()/50)
# plt.clim(-10,10)
# plt.colorbar()
# plt.show()

print("fields")
print((datad**2).mean())
print((datar**2).mean())
print(((datad-datar)**2).mean())

print(datad.max())

lim = 0.03
fig = plt.figure()
plt.imshow(datad,aspect='auto');plt.clim(-lim, lim)
plt.title('Campo direto em t=100')
fig.savefig('fieldd')

fig = plt.figure()
plt.imshow(datar,aspect='auto');plt.clim(-lim, lim)
plt.title('Campo reconstruido em t=100')
fig.savefig('fieldr')

fig = plt.figure()
plt.imshow(datad-datar,aspect='auto');plt.clim(-lim, lim)
plt.title('Diferença')
fig.savefig('diff')