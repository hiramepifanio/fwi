import numpy as np
import matplotlib.pyplot as plt
import random as rd

file_name = 'grad'
folderPath = "../data/"

grad = np.loadtxt(folderPath + file_name+'.data', dtype=np.float32)
grad_shape = np.loadtxt(folderPath + file_name+'.shp', dtype=int)
grad = grad.reshape(grad_shape)

gradr = np.loadtxt(folderPath + file_name+'r.data', dtype=np.float32)
gradr = gradr.reshape(grad_shape)

diff = grad - gradr
d = abs((grad - gradr)/grad)

figname = 'grads-diff-percent'
fig = plt.figure(figname, figsize=(7, 9))

plt.subplot(311)
plt.imshow(d)
plt.xlabel('$j$')
plt.ylabel('$i$')
plt.clim(0,100)
cbar = plt.colorbar(mappable = None)
cbar.set_label('Erro (%)')

plt.subplot(312)
plt.imshow(d)
plt.xlabel('$j$')
plt.ylabel('$i$')
plt.clim(0,10)
cbar = plt.colorbar(mappable = None)
cbar.set_label('Erro (%)')

plt.subplot(313)
plt.imshow(d)
plt.xlabel('$j$')
plt.ylabel('$i$')
plt.clim(0,1)
cbar = plt.colorbar(mappable = None)
cbar.set_label('Erro (%)')


fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)

plt.show()
fig.savefig(figname)

# d = 100*diff/grad
# erro = np.mean(d)
# print(erro)

# erro = 0
# for i in range(2, grad_shape[0] - 2):
#     for j in range(2, grad_shape[1] - 2):
#         erro += 100*((((grad[i,j] - gradr[i,j])/grad[i,j])**2)**0.5)
# erro /= (grad_shape[0]*grad_shape[1])
# print(erro)

# def p(a):
#     return format(a,'.1E')

# for r in range(20):
#     i = int(rd.random()*(grad_shape[0]))
#     j = int(rd.random()*(grad_shape[1]))
#     a, b = grad[i, j], gradr[i, j]
#     print(i, j, p(a), p(b), p(a - b), p(abs((a-b)/a)), p(100*abs((a-b)/a)))

# d = abs((grad - gradr)/grad)
# for i in range(grad_shape[0]):
#     print(i, d[i,:].max())