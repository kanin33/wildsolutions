import time
import numpy as np
import matplotlib.pyplot as plt
from streamfunction import plot_streamlines, plot_curl
from user_cmap import MidpointNormalize, add_colorbar

def movie_particles(v, x1_min, x1_max,
                    x2_min, x2_max,
                    t_min, t_max,
                    n, frames,
                    mu_0, mu_1, trace=1,
                    output="movie_particles.gif",
                    facecolor="k",
                    cmap="viridis"):
    """
    Make a movie visualizing the vector field v(t, x1, x2) in time.

    :param v:      two dimensional vector field of (t, x1, x2)
    :param x1_min: minimum x1 value
    :param x1_max: maximum x1 value
    :param t_min:  minimum x2 value
    :param t_max:  maximum x2 value
    :param n:      number of particles
    :param frames: number of frames
    :param mu_0:   x2_min = t*mu_0
    :param mu_1:   x2_max = t*mu_1
    :param trace:  number of points backwards in time to be shown
    """
    import os
    # make particles
    px, py    = np.mgrid[x1_min:x1_max:n*1j, x2_min:x2_max:n*1j]
    particles = [(px, py)]*trace
    dt = float(t_max-t_min)/(frames)

    plt.clf()
    for i in range(frames):
        print "t=%g" % (dt*i + t_min)
        U, V = runge_kutta_4(v, dt, py, px, t_min+i*dt)
        px  += dt*U
        py  += dt*V
        for j in range(trace-1):
            px_j, py_j   = particles[j+1]
            particles[j] = (np.copy(px_j), np.copy(py_j))
        particles[trace-1] = (px, py)

        colors = np.sqrt(U**2 + V**2)

        plt.gca().set_aspect('equal', adjustable='box')
        plt.axis([x2_min, x2_max, x1_min, x1_max])
        for j in range(trace):
            px_j, py_j = particles[j]
            alpha=(0.01 + 0.9*(j+1)/float(trace))
            plt.scatter(py_j, px_j, c=colors, s=0.05, alpha=alpha, cmap=cmap) # marker=(5, 1, 10),
        ax = plt.gca()
        import matplotlib as mpl
        if int(mpl.__version__[0]) <= 1:
            ax.set_axis_bgcolor("k")
        else:
            ax.set_facecolor(facecolor)
        plt.xticks([])
        plt.yticks([])
        plt.savefig('frame_p_%04d.png' % (i+1), dpi=300)
        plt.clf()
        plt.hold("off")

    files  = 'frame_p_*.png'
    os.system("convert -delay 10 -layers OptimizePlus %s %s" % (files, output))

def runge_kutta_4(v, dt, px, py, t_min):
    UK1, VK1 = v(np.float64(t_min), py, px)
    UK2, VK2 = v(np.float64(t_min + dt/2.), py + dt*VK1/2, px + dt*UK1/2)
    UK3, VK3 = v(np.float64(t_min + dt/2.), py + dt*VK2/2, px + dt*UK2/2)
    UK4, VK4 = v(np.float64(t_min + dt), py + dt*VK3, px + dt*UK3)
    return UK4, VK4

def plot_v_x1_x2(v, x2_min, x2_max, x1_min, x1_max, t, n=50,
                 stream=False, contour=False, heatmap=False, quiver=False,
                 stream_plt=False, curl=False,
                 timing=True, filename="testpicture.png",
                 show=True,
                 dpi=None,
                 num_levels=30):
    """
    Plot the vector field given by v(t, x1, x2).

    Choose a plot methond among
    stream  - streamplot
    contour - contourplot of size of vectors
    heatmap - heatmap of size of vectors
    quiver  - quiver plot
    stream_plt - matplotlib streamplot
    curl    - plot the curl
    """
    if stream:
        plot_streamlines(lambda x, y: v(t, x, y),
                         x1_min, x1_max,
                         x2_min, x2_max,
                         n=500,
                         num_levels=num_levels)
        plt.savefig(filename, format='png', dpi=dpi)
        plt.hold(False)
        if show:
            plt.show()
        return
    if curl:
        plot_curl(lambda x, y: v(t, x, y),
                  x1_min, x1_max,
                  x2_min, x2_max,
                  n=n)
        plt.savefig(filename, format='png', dpi=dpi)
        plt.hold(False)
        return

    X, Y = np.mgrid[x1_min:x1_max:n*1j, x2_min:x2_max:n*1j]
    x    = np.linspace(x1_min, x1_max, n)
    y    = np.linspace(x2_min, x2_max, n)

    t0   = time.time()
    U, V = v(t, X, Y)
    t1   = time.time()

    if timing:
        print "Time used to compute %g x %g mesh: %f" % (n, n, (t1 - t0))

    Z    = np.sqrt(U**2 + V**2)

    # Choose a color pallette. The lowermost will be used
    cmap = plt.cm.autumn
    cmap = plt.cm.Purples
    cmap = plt.cm.BuPu
    cmap = plt.cm.Greys
    cmap = plt.cm.Blues
    plt.grid(False)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis([x2_min, x2_max, x1_min, x1_max])

    if stream_plt:
        plt.streamplot(y, x, U, V, color=U, cmap=cmap,
                       density=6,
                       linewidth=0.5,
                       arrowstyle="-",
                       minlength=0.001)
    if contour:
        cs = plt.contourf(y, x, Z, cmap=cmap)
        plt.colorbar(cs, shrink=0.8)
    if quiver:
        arrows_x = int(n)/50
        arrows_y = int(n)/50
        plt.quiver(Y[::arrows_x,::arrows_y], X[::arrows_x,::arrows_y],
                   U[::arrows_x,::arrows_y], V[::arrows_x,::arrows_y])

    plt.savefig(filename, format='png', dpi=dpi)
    if show:
        plt.show()

def plot_scalar_field(v, x1_min, x1_max, x2_min, x2_max, n,
                      cmap=plt.cm.BuPu,
                      zero_value=0, f_min=None, f_max=None,
                      title="",
                      subplot=111,
                      shrink=1,
                      ticks=None):
    """
    Plot scalar field of x and y given by v
    """
    X, Y = np.mgrid[x1_min:x1_max:n*1j, x2_min:x2_max:n*1j]
    f = v(X, Y)
    #cmap = plt.cm.autumn
    #cmap = plt.cm.Purples
    #cmap = plt.cm.BuPu
    #cmap = plt.cm.Greys
    #cmap = plt.cm.Blues
    fig = plt.gcf()
    ax = fig.add_subplot(subplot)
    plt.grid(False)
    plt.title(title)
    plt.xticks([], [])
    plt.yticks([], [])
    if f_min is None:
        f_min = np.min(f)
    if f_max is None:
        f_max = np.max(f)
    im = ax.imshow(f, cmap=cmap, clim=(f_min, f_max),
                   norm=MidpointNormalize(midpoint=zero_value, vmin=f_min, vmax=f_min))
    plt.gca().set_aspect((x1_max-x1_min)/(x2_max-x2_min))
    plt.colorbar(im, shrink=shrink, ticks=ticks)
    #plt.show()

def movie_scalar_field(v, x1_min, x1_max, x2_min, x2_max, t_min, t_max, n,
                       output="scalar_movie.gif",
                       frames=30,
                       cmap=plt.cm.BuPu,
                       zero_value=0, f_min=None, f_max=None,
                       title="",
                       subplot=111,
                       shrink=1,
                       ticks=None):

    import os
    dt = (t_max - t_min)/(frames-1)
    for i in range(frames):
        X, Y = np.mgrid[x1_min:x1_max:n*1j, x2_min:x2_max:n*1j]
        f = v(t_min + dt*i, X, Y)
        #cmap = plt.cm.autumn
        #cmap = plt.cm.Purples
        #cmap = plt.cm.BuPu
        #cmap = plt.cm.Greys
        #cmap = plt.cm.Blues
        fig = plt.gcf()
        ax = fig.add_subplot(subplot)
        plt.grid(False)
        plt.title(title + " t=%f" % (t_min + dt*i))
        plt.xticks([], [])
        plt.yticks([], [])
        if f_min is None:
            f_min = np.min(f)
        if f_max is None:
            f_max = np.max(f)
        im = ax.imshow(f, cmap=cmap, clim=(f_min, f_max),
                       norm=MidpointNormalize(midpoint=zero_value, vmin=f_min, vmax=f_min))
        plt.gca().set_aspect((x1_max-x1_min)/(x2_max-x2_min))
        plt.colorbar(im, shrink=shrink, ticks=ticks)
        plt.savefig('frame_p_%04d.png' % (i+1), dpi=300)
        plt.clf()
        plt.hold("off")

    files  = 'frame_p_*.png'
    os.system("convert -delay 10 -layers OptimizePlus %s %s" % (files, output))
