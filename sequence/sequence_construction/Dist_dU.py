from Grid import Grid

class Dist_dU:
    """
    Class for finding the distance from a point to tbe boundary of dU.
    The distance values should be precomputed and stored in a file
    """

    def __init__(self, dist_filename):
        """
        Read the given file containing distance values

        File should be on the form
        n xmin xmax
        llllkkkkjjjjiiii value
        """
        with open(dist_filename, 'r') as infile:
            n, xmin, xmax = infile.readline().split()
            self.dist = Grid(0, [int(n)]*4)
            self.xmin = float(xmin)
            self.xmax = float(xmax)
            self.dx   = (self.xmax - self.xmin)/float(n)
            for line in infile:
                index, value = line.split()
                i = int(index[-3:])
                j = int(index[-6:-3])
                k = int(index[-9:-6])
                l = int(index[:-9])
                self.dist[(i, j, k, l)] = float(value)

    def __call__(self, point):
        """
        Use linear interpolation between surrounding values
        :param point: [x, y, z, q]
        :return: distance to boundary
        """
        x = point[0]; y = point[1]; z = point[2]; q = point[3]
        x_ind = (x - self.xmin)/self.dx
        y_ind = (y - self.xmin)/self.dx
        z_ind = (z - self.xmin)/self.dx
        q_ind = (q - self.xmin)/self.dx
        i = int(x_ind)
        j = int(y_ind)
        k = int(z_ind)
        l = int(q_ind)
        x_rest = x_ind - i
        y_rest = y_ind - j
        z_rest = z_ind - k
        q_rest = q_ind - l
        dist_ip1 = abs(self.dist[(i+1, j, k, l)])
        dist_jp1 = abs(self.dist[(i, j+1, k, l)])
        dist_kp1 = abs(self.dist[(i, j, k+1, l)])
        dist_lp1 = abs(self.dist[(i, j, k, l+1)])
        dist     = abs(self.dist[(i, j, k, l)])
        x_interpol = dist_ip1*x_rest + dist*(1 - x_rest)
        y_interpol = dist_jp1*y_rest + dist*(1 - y_rest)
        z_interpol = dist_kp1*z_rest + dist*(1 - z_rest)
        q_interpol = dist_lp1*q_rest + dist*(1 - q_rest)

        return (x_interpol + y_interpol + z_interpol + q_interpol)/4.

if __name__ == '__main__':
    d = Dist_dU("lol")
    print d([-0.000001, 0, 0, 0])
    print d([0.1, 0, 0, 0])
    print d([-3, 3, 4, 5])

