import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import lil_matrix, kronsum
from scipy.sparse.linalg import spsolve

def plot_streamlines(v, xmin, xmax, ymin, ymax, n, num_levels=30,
                     streamplot=False):
    """
    Plot streamlines of vector field v = [U, V]

    :param v: vector field of x and y
    :param xmin:
    :param xmax:
    :param ymin:
    :param ymax:
    :param n: mesh size
    """

    T1 = lil_matrix((n+2,n+2))
    T1[0,0] = 1
    T1[n+1, n+1] = 1
    for i in range(1,n+1):
        T1[i,i] = 2
    for i in range(n+1):
        T1[i+1,i] = -1
        T1[i, i+1] = -1

    T2 = kronsum(T1, T1).tolil()

    # Set one point to a constant
    T2[0,0] = 1
    T2[0,1] = 0
    T2[0,n+2] = 0
    T2 = T2.tocsr()

    h = (xmin-xmax)/float(n+2)

    X, Y = np.mgrid[xmin:xmax+h:(n+3)*1j, ymin:ymax+h:(n+3)*1j]

    U, V = v(X, Y)

    b = (-(V[1:,:-1] - V[1:,1:]) + (U[:-1,1:] - U[1:,1:]))*h

    # Neumann boundary condition
    b[:,0]   += +h*V[:-1,0]
    b[:,n+1] += -h*V[1:,n+1]
    b[0,:]   += -h*U[0,:-1]
    b[n+1,:] += +h*U[n+1,1:]

    # Set one point to a constant
    b[0,0] = 0.5

    b = b.flatten()

    u   = spsolve(T2, b)
    phi = u.reshape((n+2, n+2))
    x   = np.linspace(xmin, xmax, n+2)
    y   = np.linspace(ymin, ymax, n+2)
    phimin = np.min(phi)
    phimax = np.max(phi)
    diff = phimax - phimin
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis([ymin, ymax, xmin, xmax])
    levels = np.linspace(phimin + diff/10., phimax - diff/10.*2, num_levels)
    plt.contour(y, x, phi, levels=levels, linewidths=1, cmap=plt.cm.Greys)# cmap=plt.cm.Blues)
    plt.xticks([], [])
    plt.yticks([], [])

    if streamplot:
        x = np.linspace(xmin, xmax, n+3)
        y = np.linspace(ymin, ymax, n+3)
        plt.streamplot(x, y, U, V)
    #plt.show()

def plot_curl(v, xmin, xmax, ymin, ymax, n):
    h1 = (xmin-xmax)/float(n+2)
    h2 = (ymin-ymax)/float(n+2)

    X, Y = np.mgrid[xmin:xmax+h1:(n+3)*1j, ymin:ymax+h2:(n+3)*1j]

    U, V = v(X, Y)

    b = -(V[1:,:-1] - V[1:,1:])/h1 + (U[:-1,1:] - U[1:,1:])/h2
    plt.imshow(b, cmap=plt.cm.Blues, vmin=-20, vmax=20)
    plt.colorbar()
    plt.xticks([], [])
    plt.yticks([], [])
    plt.show()

def plot_circular_streamlines(n):
    # Test for streamfunction with V = [2y, -2x]
    v = lambda x, y: [2*y, -2*x]
    plot_streamlines(v, -1, 1, -1, 1, n)

def plot_circular_curl(n):
    v = lambda x, y: [2*y, -2*x]
    plot_curl(v, -1, 1, -1, 1, n)

if __name__ == '__main__':
    print ""
    plot_circular_curl(10)
    #plot_circular_streamlines(1)
    #plot_ex_5_6_streamlines(200)
