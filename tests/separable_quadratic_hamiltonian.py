# Module that implements functions for the separable quadratic normal form Hamiltonian with index one saddle


## Plotting quadratic Hamiltonian (saddle projection)

import numpy as np

from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt


def normal_form_potential(q1, q2, LAMBDA = 1, OMEGA2 = 1):
    """ Potential energy function for the uncoupled normal form Hamiltonian
    """

    normal_form_pe = -0.5*LAMBDA*q1**2 + 0.5*OMEGA2*q2**2 

    return normal_form_pe


# Plotting the potential energy function
x = np.arange(-1.0, 1.0, 0.01)
y = np.arange(-1.0, 1.0, 0.01)
X,Y = np.meshgrid(x, y) # grid of point
Z = normal_form_potential(X, Y) # evaluation of the function on the grid

fig = plt.figure()
ax = fig.gca()
im = ax.imshow(Z, cmap = cm.RdBu) # drawing the function
# ax.set_xlim([-1,1])
# ax.set_ylim(([-1,1])

# adding the Contour lines with labels
cset = contour(X,Y,Z, np.arange(0, 0.4, 40), linewidths = 2, cmap = cm.Set2)
clabel(cset, inline = True,fmt='%1.2f', fontsize=10)
colorbar(im, shrink = 0.7, pad = 0.02,
        drawedges = True) 
show()


fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride = 2, cstride = 2, 
                       cmap = cm.RdBu, linewidth = 0, antialiased=True)

# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()





