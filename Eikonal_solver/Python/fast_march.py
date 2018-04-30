import numpy as np
from initialize import initialize_4D
from heapq import heapify, heappush, heappop
from Grid import Grid

def fast_march(grid, xmin, xmax, ymin, ymax, zmin, zmax, qmin, qmax, boundary, F):
    # Initialize
    grid, narrow_band, not_known = initialize_4D(grid, xmin, xmax,
                                                       ymin, ymax,
                                                       zmin, zmax,
                                                       qmin, qmax,
                                                       boundary)
    print "Done with initialize"

    # March forward
    grid = march_forward(grid, xmin, xmax, ymin, ymax, zmin, zmax, qmin, qmax,
                         narrow_band, not_known, F)
    print "Done with fast march"

    return grid

def march_forward(grid, xmin, xmax, ymin, ymax, zmin, zmax, qmin, qmax,
                  narrow_band, not_known, F):
    n1, n2, n3, n4 = grid.shape
    dx = (xmax - xmin) / float(n1 - 1)
    dy = (ymax - ymin) / float(n2 - 1)
    dz = (zmax - zmin) / float(n3 - 1)
    dq = (qmax - qmin) / float(n4 - 1)
    heapify(narrow_band)

    while not_known > 0:
        # step 1
        A = heappop(narrow_band)

        # step 2
        i, j, k, l = A[1], A[2], A[3], A[4]
        grid[i, j, k, l] = -A[0]

        # step 3
        gip1 = grid[i+1, j, k, l]
        # If Trial or Far, update value
        if gip1 > 0:
            grid[i+1, j, k, l] = advance(grid, i+1, j, k, l, dx, dy, dz, dq,
                                      F, xmin, ymin, zmin, qmin)
            # If far, add to narrow band
            if gip1 == np.inf:
                heappush(narrow_band, (grid[i+1, j, k, l], i+1, j, k, l))

        gip1 = grid[i-1, j, k, l]
        if gip1 > 0:
            grid[i-1, j, k, l] = advance(grid, i-1, j, k, l, dx, dy, dz, dq,
                                      F, xmin, ymin, zmin, qmin)
            if gip1 == np.inf:
                heappush(narrow_band, (grid[i-1, j, k, l], i-1, j, k, l))

        gip1 = grid[i, j+1, k, l]
        if gip1 > 0:
            grid[i, j+1, k, l] = advance(grid, i, j+1, k, l, dx, dy, dz, dq,
                                      F, xmin, ymin, zmin, qmin)
            if gip1 == np.inf:
                heappush(narrow_band, (grid[i, j+1, k, l], i, j+1, k, l))

        gip1 = grid[i, j-1, k, l]
        if gip1 > 0:
            grid[i, j-1, k, l] = advance(grid, i, j-1, k, l, dx, dy, dz, dq,
                                      F, xmin, ymin, zmin, qmin)
            if gip1 == np.inf:
                heappush(narrow_band, (grid[i, j-1, k, l], i, j-1, k, l))

        gip1 = grid[i, j, k+1, l]
        if gip1 > 0:
            grid[i, j, k+1, l] = advance(grid, i, j, k+1, l, dx, dy, dz, dq,
                                      F, xmin, ymin, zmin, qmin)
            if gip1 == np.inf:
                heappush(narrow_band, (grid[i, j, k+1, l], i, j, k+1, l))

        gip1 = grid[i, j, k-1, l]
        if gip1 > 0:
            grid[i, j, k-1, l] = advance(grid, i, j, k-1, l, dx, dy, dz, dq,
                                      F, xmin, ymin, zmin, qmin)
            if gip1 == np.inf:
                heappush(narrow_band, (grid[i, j, k-1, l], i, j, k-1, l))

        gip1 = grid[i, j, k, l+1]
        if gip1 > 0:
            grid[i, j, k, l+1] = advance(grid, i, j, k, l+1, dx, dy, dz, dq,
                                         F, xmin, ymin, zmin, qmin)
            if gip1 == np.inf:
                heappush(narrow_band, (grid[i, j, k, l+1], i, j, k, l+1))

        gip1 = grid[i, j, k, l-1]
        if gip1 > 0:
            grid[i, j, k, l-1] = advance(grid, i, j, k, l-1, dx, dy, dz, dq,
                                         F, xmin, ymin, zmin, qmin)
            if gip1 == np.inf:
                heappush(narrow_band, (grid[i, j, k, l-1], i, j, k, l-1))
        not_known -= 1

    return grid

def advance(grid, i, j, k, l, dx, dy, dz, dq, F, xmin, ymin, zmin, qmin):
    F_ij = F(i*dx + xmin, j*dy + ymin, k*dz + zmin, l*dq + qmin)
    T_H = min(abs(grid[i+1, j, k, l]), abs(grid[i-1, j, k, l]))
    T_V = min(abs(grid[i, j+1, k, l]), abs(grid[i, j-1, k, l]))
    T_Z = min(abs(grid[i, j, k+1, l]), abs(grid[i, j, k-1, l]))
    T_Q = min(abs(grid[i, j, k, l+1]), abs(grid[i, j, k, l-1]))
    T = [T_H, T_V, T_Z, T_Q]
    T = [T_j for T_j in T if T_j != np.inf]
    T.sort()
    j = len(T)
    while sum([(T[j-1] - T[i])**2 for i in range(j-1)]) >= dx**2 / float(F_ij**2):
        j = j - 1
    return 1./j*sum([T[i] for i in range(j)])\
         + 1./j*np.sqrt((sum([T[i] for i in range(j)]))**2
         - j*(sum([T[i]**2 for i in range(j)]) - dx**2 / float(F_ij**2)))

def exact_square(n, xmin, xmax, r):
    grid = -np.zeros((n, n, n, n))
    h = (xmax - xmin)/float(n-1)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    x = abs(i*h + xmin)
                    y = abs(j*h + xmin)
                    z = abs(k*h + xmin)
                    q = abs(l*h + xmin)
                    if x > 0.9 or y > 0.9 or z > 0.9 or q > 0.9:
                        grid[i,j,k,l] = -0
                    else:
                        grid[i,j,k,l] = -(0.9 - max(x, y, z, q))
    return grid

def exact_circle(n, xmin, xmax, r):
    grid = -np.zeros((n, n, n, n))
    h = (xmax - xmin)/float(n-1)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    x = abs(i*h + xmin)
                    y = abs(j*h + xmin)
                    z = abs(k*h + xmin)
                    q = abs(l*h + xmin)
                    if x**2 + y**2 + z**2 + q**2 > r**2:
                        grid[i,j,k,l] = -0
                    else:
                        grid[i,j,k,l] = -(0.9 - np.sqrt(x**2 + y**2 + z**2 + q**2))
    return grid

def compute_dist_dU(C=1, m=10):
    # boundary used in the thesis
    n1, n2, n3, n4 = [m]*4
    grid = -np.zeros((n1, n2, n3, n4))
    xmin, xmax = -np.sqrt(C), np.sqrt(C)
    ymin, ymax = -np.sqrt(C), np.sqrt(C)
    zmin, zmax = -np.sqrt(C), np.sqrt(C)
    qmin, qmax = -np.sqrt(C), np.sqrt(C)
    def boundary(x, y, z, q):
        # compute both eigenvalues
        lambda_1 = 1/2. * (-(x**2 + y**2) + C) + 1/2. \
                  * np.sqrt(x**4 + 2*x**2*y**2 - 4*x**2*z - 8*x*y*q
                            + y**4 + 4*y**2*z + 4*z**2 + 4*q**2)
        lambda_2 = 1/2. * (-(x**2 + y**2) + C) - 1/2. \
                  * np.sqrt(x**4 + 2*x**2*y**2 - 4*x**2*z - 8*x*y*q
                            + y**4 + 4*y**2*z + 4*z**2 + 4*q**2)
        if lambda_1 > 0 and lambda_2 > 0:
            return -lambda_1*lambda_2
        return 0
    F = lambda x, y, z, q: 1
    grid = fast_march(grid, xmin, xmax, ymin, ymax, zmin, zmax, qmin, qmax,
                      boundary, F)
    return grid

def write_to_file(grid, xmin, xmax, filename="python_dist.txt"):
    n = grid.shape[0]
    with open(filename, "w") as outfile:
        outfile.write("%g %g %g \n" % (n, xmin, xmax))
        for point in grid.elements:
            outfile.write("%03d%03d%03d%03d " % (point))
            outfile.write("%g\n" % grid[point])
    return


if __name__ == '__main__':
    solution_grid = Grid(standard_value=0, shape=[30] * 4)

    def boundary_function(x, y, z, q):
        return x ** 2 + y ** 2 + z ** 2 + q ** 2 - 1

    fast_march(grid=solution_grid,
               xmin=-2, xmax=2,
               ymin=-2, ymax=2,
               zmin=-2, zmax=2,
               qmin=-2, qmax=2,
               boundary=boundary_function,
               F=lambda x, y, z, q: 1)
    print solution_grid[15, 15, 15, 15]
