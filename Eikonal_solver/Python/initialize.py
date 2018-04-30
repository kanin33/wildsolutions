import numpy as np

def initialize_4D(grid, xmin, xmax, ymin, ymax, zmin, zmax, qmin, qmax, boundary):
    boundary_bool = lambda x, y, z, q: boundary(x, y, z, q) < 0
    n1, n2, n3, n4 = grid.shape
    dx   = (xmax - xmin) / float(n1 - 1)
    dy   = (ymax - ymin) / float(n2 - 1)
    dz   = (zmax - zmin) / float(n3 - 1)
    dq   = (qmax - qmin) / float(n4 - 1)
    # initialize inner points
    not_known = 0
    inner_points = []
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                for l in range(n4):
                    if boundary_bool(i*dx + xmin, j*dy + ymin,
                                     k*dz + zmin, l*dq + qmin):
                        grid[i,j,k,l] = np.inf
                        not_known += 1
                        inner_points.append([i, j, k, l])
    # initialize narrow band
    # assume boundary is inside domain
    # store the narrow points as list of tuples (value, i, j)
    narrow_band = []
    for inner_point in inner_points:
        # s1 = d(p, pi+1)
        # s2 = d(p, pj+1)
        # s3 = d(p, pk+1)
        # s4 = d(p, pi-1)
        # s5 = d(p, pj-1)
        # s6 = d(p, pk-1)
        i, j, k, l = inner_point
        s1, s2, s3, s4, s5, s6, s7, s8 = [np.inf]*8
        if grid[i+1, j, k, l] == 0:
            p0 = np.array([xmin + i*dx, ymin + j*dy, zmin + k*dz, qmin + l*dq])
            p1 = np.array([xmin + (i+1)*dx, ymin + j*dy, zmin + k*dz, qmin + l*dq])
            s1 = dist_p0_boundary(p0, p1, boundary)
        if grid[i-1, j, k, l] == 0:
            p0 = np.array([xmin + i*dx, ymin + j*dy, zmin + k*dz, qmin + l*dq])
            p1 = np.array([xmin + (i-1)*dx, ymin + j*dy, zmin + k*dz, qmin + l*dq])
            s5 = dist_p0_boundary(p0, p1, boundary)
        if grid[i, j+1, k, l] == 0:
            p0 = np.array([xmin + i*dx, ymin + j*dy, zmin + k*dz, qmin + l*dq])
            p1 = np.array([xmin + i*dx, ymin + (j+1)*dy, zmin + k*dz, qmin + l*dq])
            s2 = dist_p0_boundary(p0, p1, boundary)
        if grid[i, j-1, k, l] == 0:
            p0 = np.array([xmin + i*dx, ymin + j*dy, zmin + k*dz, qmin + l*dq])
            p1 = np.array([xmin + i*dx, ymin + (j-1)*dy, zmin + k*dz, qmin + l*dq])
            s6 = dist_p0_boundary(p0, p1, boundary)
        if grid[i, j, k+1, l] == 0:
            p0 = np.array([xmin + i*dx, ymin + j*dy, zmin + k*dz, qmin + l*dq])
            p1 = np.array([xmin + i*dx, ymin + j*dy, zmin + (k+1)*dz, qmin + l*dq])
            s3 = dist_p0_boundary(p0, p1, boundary)
        if grid[i, j, k-1, l] == 0:
            p0 = np.array([xmin + i*dx, ymin + j*dy, zmin + k*dz, qmin + l*dq])
            p1 = np.array([xmin + i*dx, ymin + j*dy, zmin + (k-1)*dz, qmin + l*dq])
            s7 = dist_p0_boundary(p0, p1, boundary)
        if grid[i, j, k, l+1] == 0:
            p0 = np.array([xmin + i*dx, ymin + j*dy, zmin + k*dz, qmin + l*dq])
            p1 = np.array([xmin + i*dx, ymin + j*dy, zmin + k*dz, qmin + (l+1)*dq])
            s4 = dist_p0_boundary(p0, p1, boundary)
        if grid[i, j, k, l-1] == 0:
            p0 = np.array([xmin + i*dx, ymin + j*dy, zmin + k*dz, qmin + l*dq])
            p1 = np.array([xmin + i*dx, ymin + j*dy, zmin + k*dz, qmin + (l-1)*dq])
            s8 = dist_p0_boundary(p0, p1, boundary)
        t1 = dx*min(s1, s5)
        t2 = dy*min(s2, s6)
        t3 = dz*min(s3, s7)
        t4 = dq*min(s4, s8)
        grid[i,j,k,l] = np.sqrt(np.float64(1)/float(1/t2**2 + 1/t1**2 + 1/t3**2 + 1/t4**2))
        if grid[i,j,k,l] != np.inf:
            narrow_band.append((grid[i,j,k,l], i, j, k, l))

    return grid, narrow_band, not_known

def dist_p0_boundary(p0, p1, f):
    t = f(*p0)/float(f(*p1) - f(*p0))
    return length(t*(p0 - p1))

def length(vector):
    return np.sqrt(np.sum(vector**2))
