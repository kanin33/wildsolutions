from sequence_functions import *
import os

class Sequence:
    """
    Class for representing the sequence of functions.
    """
    def __init__(self, C, P, dist_obj, v_tilde, u_tilde, dirname="test_dir"):
        """
        :param C:        constant
        :param P:        domain, [t_min, t_max, mu_0, mu_1, x1_min, x1_max]
        :param dist_obj: object that takes a point (t, x1, x2) and gives distance
        :param v_tilde:  constant [v1, v2]
        :param u_tilde:  constant [[u1, u2], [u2, -u1]]
        """
        self.C = C

        self.P       = P
        self.t_min   = P[0]
        self.t_max   = P[1]
        self.mu_0    = P[2]
        self.mu_1    = P[3]
        self.x_1_min = P[4]
        self.x_1_max = P[5]
        self.P_size  = self.size_P()

        self.dist_obj = dist_obj

        self.v_tilde = v_tilde
        self.u_tilde = u_tilde

        self.r_values     = [] # list with all r-valies
        self.num_points_r = [] # number of points with radius i
        self.point_array  = np.zeros((0, 11 + 1 + 5 + 1))


        v_tilde  = self.v_tilde
        P_size   = self.P_size

        self.integral = (v_tilde[0]**2 + v_tilde[1]**2)*P_size
        self.integral_prev = 0
        self.integral_last_step = 0 # \int (z_k+1 - z_k)^2 dx

        self.dirname = dirname # name of directory to store the computed coefficients

        # make directory for storing files with coefficients
        if not os.path.exists(dirname):
            os.makedirs(dirname)

    def __call__(self, t, x1, x2):
        """
        Divide the grid into squares around each cylinder points, and
        compute the values in each points
        """
        point_array = self.point_array

        # NB: When making particle trace movies the below method must be used
        # Use if True instead
        if x1.shape == () or x2.shape == ():
        #if True:
            # if x1 or x2 are not arrays, use another call method
            return np.sum(v_p_j(t, x1, x2, *point_array.T), -1)

        f = np.array([[0*x1, 0*x1], [0*x1, 0*x1]]) # empty return array

        if len(x1.shape) == 3: # array in both x1, x2 and t direction
            min_x1 = x1[0,0,0]; min_x2 = x2[0,0,0]
            dx1 = x1[0,1,0] - min_x1
            dx2 = x2[0,0,1] - min_x2
        else:
            min_x1 = x1[0,0]; min_x2 = x2[0,0]
            dx1 = x1[1,0] - min_x1
            dx2 = x2[0,1] - min_x2

        r_values = self.r_values

        p = 0 # counter for points_array

        # iterate through all r_values
        x2_len = x1.shape[-1]
        x1_len = x1.shape[-2]
        for j in range(len(r_values[0:])):
            r = r_values[j]
            nx1 = int(r/dx1)+1
            nx2 = int(r/dx2)+1
            # for all points with radius r
            for k in range(self.num_points_r[j]):
                x, y = point_array[p][2:4]
                mid_index_x1 = int((x - min_x1)/dx1)+1
                mid_index_x2 = int((y - min_x2)/dx2)+1
                low_index_x1 = mid_index_x1 - nx1
                low_index_x2 = mid_index_x2 - nx2
                high_index_x1 = mid_index_x1 + nx1 + 1
                high_index_x2 = mid_index_x2 + nx2 + 1
                if low_index_x1 < 0:
                    low_index_x1 = 0
                if low_index_x2 < 0:
                    low_index_x2 = 0
                if high_index_x1 < 0:
                    high_index_x1 = 0
                if high_index_x2 < 0:
                    high_index_x2 = 0
                if low_index_x1 > x1_len-1:
                    low_index_x1 = x1_len-1
                if low_index_x2 > x2_len-1:
                    low_index_x2 = x2_len-1
                if high_index_x1 > x1_len-1:
                    high_index_x1 = x1_len-1
                if high_index_x2 > x2_len-1:
                    high_index_x2 = x2_len-1
                f[:,:,low_index_x1:high_index_x1, low_index_x2:high_index_x2] += \
          v_p_j(t, x1[low_index_x1:high_index_x1, low_index_x2:high_index_x2],
                   x2[low_index_x1:high_index_x1, low_index_x2:high_index_x2],
                              *point_array[p])[:,:,:,:,0]
                p += 1
        return f

    def add_points(self, points, r, n_trials, compute_integrals=True):
        """
        Add new cylinder points to the function

        :param points:   Array of points on the form [t, x1, x2]
        :param r:        Corresponding radius
        :param n_trials: Number of iterations used
        """
        p_a_size    = self.point_array.shape[0]
        point_array = np.zeros((p_a_size + len(points), 11 + 1 + 5 + 1))
        point_array[0:p_a_size, :] = self.point_array
        self.r_values.append(r)
        self.num_points_r.append(len(points))

        for i in range(self.num_points_r[-1]):
            point_array[p_a_size+i, 0] = r
            point = points[i]
            t, x, y = point
            point_array[p_a_size+i, 1] = t
            point_array[p_a_size+i, 2] = x
            point_array[p_a_size+i, 3] = y
            v_star, u_star = self.__call__(*point)
            v = v_star + self.v_tilde
            u = u_star + self.u_tilde[0,:]
            a, b, alpha_j  = find_a_b(v, u, self.C)
            eta, lambda_, p = compute_eta_p(v, u, self.C)
            N = compute_N(eta, k=len(self.r_values))
            point_array[p_a_size+i, 4] = p[0]
            point_array[p_a_size+i, 5] = p[1]
            point_array[p_a_size+i, 6] = p[2]
            point_array[p_a_size+i, 7] = p[3]
            point_array[p_a_size+i, 8] = N
            point_array[p_a_size+i, 9] = eta[0]
            point_array[p_a_size+i, 10] = eta[1]
            point_array[p_a_size+i, 11] = eta[2]
            point_array[p_a_size+i, 12] = a[0]
            point_array[p_a_size+i, 13] = a[1]
            point_array[p_a_size+i, 14] = b[0]
            point_array[p_a_size+i, 15] = b[1]
            point_array[p_a_size+i, 16] = lambda_
            point_array[p_a_size+i, 17] = 2*r # This can be changed

        self.point_array = point_array

        self.integral_prev = self.integral
        if compute_integrals:
            self.integral      = self.compute_integral_absolute()
            print "integral: ", self.integral
            print "number of points: ",
            print len(self.point_array)
            self.integral_last_step = self.compute_integral_absolute_last_step()
        else:
            self.integral = 0
            self.integral_last_step = 0

    def remove_last_points(self):
        """
        Remove the points added in last iteration from point_array
        """
        p_a_size     = self.point_array.shape[0]
        new_p_a_size = p_a_size - self.num_points_r[-1]
        point_array  = np.zeros((new_p_a_size, 11 + 1 + 5 + 1))
        point_array[:, :] = self.point_array[0:new_p_a_size, :]
        self.point_array  = point_array
        del self.num_points_r[-1]
        del self.r_values[-1]

    def compute_integral_absolute(self):
        """
        Compute int (v_tilde + v_k)^2 dx
        """
        t1 = time.time()
        v_tilde  = self.v_tilde
        P_size   = self.P_size
        integral = (v_tilde[0]**2 + v_tilde[1]**2)*P_size
        integral += integrate_matlab(v_tilde, self.point_array)
        t2 = time.time()
        print "computed integral, num_points = %g, time: %f" % (sum(self.num_points_r), (t2 - t1))
        return integral

    def compute_integral_absolute_last_step(self):
        """
        Compute int (z_k - z_k-1)^2 dx
        """
        return integrate_matlab([0, 0], self.point_array[-self.num_points_r[-1]:])

    def integral_absolute(self):
        """
        iiint_P (v_tilde[0] + v_k[0])**2 + (v_tilde[1] + v_k[1])**2 dx1dx2dt
        """
        return self.integral

    def integral_absolute_prev(self):
        """
        Return the integral of the previous step
        """
        return self.integral_prev

    def integral_absolute_last_step(self):
        """
        Return only the integral for the points added in the last step
        """
        return self.integral_last_step

    def max_grad(self, t, x1, x2):
        point_array = self.point_array
        v1_grad_b = 0
        v2_grad_b = 0
        u11_grad_b = 0
        u12_grad_b = 0
        n = 0
        for i in range(len(self.r_values)):
            n_r = self.num_points_r[i]
            v1_grad_b_i, v2_grad_b_i, u11_grad_b_i, u12_grad_b_i = grad_v_p_j_b(t, x1, x2, *point_array[n:n+n_r].T)
            v1_grad_b  += np.max(v1_grad_b_i)
            v2_grad_b  += np.max(v2_grad_b_i)
            u11_grad_b += np.max(u11_grad_b_i)
            u12_grad_b += np.max(u12_grad_b_i)
            n += n_r

        return np.max([np.sqrt(v1_grad_b), np.sqrt(v2_grad_b),
                       np.sqrt(u11_grad_b), np.sqrt(u12_grad_b)])

    def error(self, t, x1, x2):
        [v1, v2], [u1, u2] = self.__call__(t, x1, x2)
        # compute the maximum norm of the matrix with entries
        # v1^2 - u1 - C/2
        # v1*v2 - u2
        # v2^2 + u1 - C/2
        v1 += self.v_tilde[0]
        v2 += self.v_tilde[1]
        u1 += self.u_tilde[0][0]
        u2 += self.u_tilde[0][1]

        a11 = np.abs(v1**2 - u1 - self.C/2)
        a12 = np.abs(v1*v2 - u2)
        a22 = np.abs(v2**2 + u1 - self.C/2)
        m1 = np.maximum(a11, a12)
        return np.abs(np.maximum(m1, a22))

    def average_error(self, n=100):
        dt = (self.t_max - self.t_min)/float(n-1)
        error = 0
        for i in range(n):
            t = np.float64(self.t_min + i*dt)
            x_2_min = self.mu_0*(t)
            x_2_max = self.mu_1*(t)
            X, Y = np.mgrid[self.x_1_min:self.x_1_max:n*1j, x_2_min:x_2_max:n*1j]
            error += np.sum(self.error(t, X, Y))
        return error/n**3

    def maximum_error(self, n=100):
        dt = (self.t_max - self.t_min)/float(n-1)
        error = 0
        for i in range(n):
            t = np.float64(self.t_min + i*dt)
            x_2_min = self.mu_0*(t)
            x_2_max = self.mu_1*(t)
            X, Y = np.mgrid[self.x_1_min:self.x_1_max:n*1j, x_2_min:x_2_max:n*1j]
            error = max(error, np.max(self.error(t, X, Y)))
        return error

    def minimum_error(self, n=100):
        dt = (self.t_max - self.t_min)/float(n-1)
        error = np.inf
        for i in range(n):
            t = np.float64(self.t_min + i*dt)
            x_2_min = self.mu_0*(t)
            x_2_max = self.mu_1*(t)
            X, Y = np.mgrid[self.x_1_min:self.x_1_max:n*1j, x_2_min:x_2_max:n*1j]
            error = min(error, np.min(self.error(t, X, Y)))
        return error

    def eigenvalues_product(self, t, x1, x2):
        """
        Compute the product of the eigenvalues of
        v tensor v - u - C/2 Id

        When the product is 0, the value lies on the boundary
        of U
        """

        [x, y], [z, q] = self.__call__(t, x1, x2)
        C = self.C
        x += self.v_tilde[0]
        y += self.v_tilde[1]
        z += self.u_tilde[0][0]
        q += self.u_tilde[0][1]
        # compute both eigenvalues
        lambda_1 = 1/2. * (-(x**2 + y**2) + C) + 1/2. \
                        * np.sqrt(x**4 + 2*x**2*y**2 - 4*x**2*z - 8*x*y*q
                                + y**4 + 4*y**2*z + 4*z**2 + 4*q**2)
        lambda_2 = 1/2. * (-(x**2 + y**2) + C) - 1/2. \
                        * np.sqrt(x**4 + 2*x**2*y**2 - 4*x**2*z - 8*x*y*q
                                + y**4 + 4*y**2*z + 4*z**2 + 4*q**2)
        r = np.zeros(lambda_1.shape)
        w = np.logical_and(lambda_1>0, lambda_2>0) # only use positive values
        r[w] = (lambda_1*lambda_2)[w]
        not_w = np.logical_not(w)
        # if either of the eigenvalues are negative, the value is outside
        r[not_w] = -np.abs(lambda_1*lambda_2)[not_w]
        return r

    def error_tilde(self):
        return self.error(np.float64(100), np.float64(0), np.float64(0))

    def eigenvalue_tilde(self):
        return self.eigenvalues_product(np.float64(100), np.float64(0), np.float64(0))

    def dist_from_dU(self, t, x1, x2):
        """
        Get distance from boundary of U
        """
        if isinstance(x1, (float, int)):
            xyzq = self.__call__(t, x1, x2)
            xyzq = [xyzq[0,0] + self.v_tilde[0], xyzq[0,1] + self.v_tilde[1],
                    xyzq[1,0] + self.u_tilde[0][0], xyzq[1,1] + self.u_tilde[0][1]]
            return self.dist_obj(xyzq)
        else:
            d = np.zeros(x1.shape)
            for i in range(x1.shape[0]):
                for j in range(x1.shape[1]):
                    d[i, j] = self.dist_from_dU(t, x1[i, j], x2[i, j])
            return d

    def size_P(self):
        """
        Compute size of P
        """
        t_min   = self.t_min
        t_max   = self.t_max
        mu_0    = self.mu_0
        mu_1    = self.mu_1
        x_1_min = self.x_1_min
        x_1_max = self.x_1_max

        x2_min_l = t_min*mu_0
        x2_max_l = t_max*mu_0
        x2_min_r = t_min*mu_1
        x2_max_r = t_max*mu_1
        x_1_len  = x_1_max - x_1_min
        return x_1_len*(x2_max_r - x2_max_l + x2_min_r - x2_min_l)*(t_max - t_min)/2.

    def v_hat_u_hat(self, t, x1, x2):
        """
        This function computes v_hat_u_hat which is an approximation of (v, u)

        Divide the grid into squares around each cylinder points, and
        compute the values in each points
        """
        if x1.shape == () or x2.shape == ():
            "if x1 or x2 are not arrays, use another call method"
            return np.sum(v_p_j_hat(t, x1, x2, *self.point_array.T), -1)

        point_array = self.point_array
        f = np.array([[0*x1, 0*x1], [0*x1, 0*x1]]) # empty return array

        dx1 = x1[1,0] - x1[0,0]
        dx2 = x2[0,1] - x2[0,0]

        t_min   = self.t_min
        t_max   = self.t_max
        mu_0    = self.mu_0
        mu_1    = self.mu_1
        x_1_min = self.x_1_min
        x_1_max = self.x_1_max
        r_values = self.r_values

        p = 0 # counter for points_array
        for r in r_values:
            h     = t_max - t_min     # height of set in t direction
            n_t   = int(h/(2*r))      # number of cylinders in t direction
            l_x_1 = x_1_max - x_1_min # length of set in x_1 direction
            n_x_1 = int(l_x_1/(2*r))  # number of cylinders in x_1 direction

            for i in xrange(n_t):
                t_i     = t_min + i*2*r
                x_2_min = mu_0*t_i
                x_2_max = mu_1*t_i
                l_x_2   = x_2_max - x_2_min
                n_x_2   = int(l_x_2/(2*r))
                # number of points on a line with length r
                nx1 = int(r/dx1)+1
                nx2 = int(r/dx2)+1
                for j in range(n_x_1):
                    for k in range(n_x_2):
                        x, y = point_array[p][2:4]
                        index_j = int((x - x1[0,0])/dx1)+1
                        index_k = int((y - x2[0,0])/dx2)+1
                        f[:,:,index_j-nx1:index_j+nx1+1, index_k-nx2:index_k+nx2+1] += \
                            v_p_j_hat(t, x1[index_j-nx1:index_j+nx1+1, index_k-nx2:index_k+nx2+1],
                                  x2[index_j-nx1:index_j+nx1+1, index_k-nx2:index_k+nx2+1],
                                  *point_array[p])[:,:,:,:,0]
                        p += 1
        return f

    def read_data_from_file(self):
        """
        Read all files with name "step_number.txt" to points_array

        File format:
        step_number
        radius
        number of points
        nu_k
        point_array-row
        ...
        point_array_row
        """
        integral_absolute = self.integral_absolute()
        print "int", integral_absolute
        current_folder = os.getcwd().split(os.sep)[-1]
        if not current_folder == self.dirname:
            os.chdir(self.dirname)
        for filename in os.listdir("."):
            if filename.startswith("step_") and filename.endswith(".txt"):
                with open(filename) as infile:
                    k = infile.readline()
                    r_k = infile.readline()
                    n_1 = infile.readline()
                    infile.readline()
                    integral_absolute = float(infile.readline())
                    k = int(k); r_k = float(r_k); n_1 = int(n_1)
                    p_a_size    = self.point_array.shape[0]
                    point_array = np.zeros((p_a_size + n_1, 11 + 1 + 5 + 1))
                    point_array[0:p_a_size, :] = self.point_array
                    self.r_values.append(r_k)
                    self.num_points_r.append(n_1)
                    i = 0
                    for line in infile:
                        values = line.split()
                        point_array[p_a_size+i, 17] = 1
                        for j in range(len(values)):
                            point_array[p_a_size+i, j] = values[j]
                        i += 1
                self.point_array = point_array
        self.integral = integral_absolute # set integral to integral in last step

    def write_step_to_file(self, nu_k):
        """
        Write values from current step to file

        File format:
        step_number
        radius
        number of points
        nu_k
        integral_absolute
        point_array-row
        ...
        point_array_row
        """
        k = len(self.r_values)            # current step
        r_k = self.r_values[-1]           # current radius
        n_1 = self.num_points_r[-1]       # number of new points
        integral_absolute = self.integral # current integral
        filename = "step_%03d.txt" % k
        current_folder = os.getcwd().split(os.sep)[-1]
        print current_folder
        if not current_folder == self.dirname:
            os.chdir(self.dirname)
        with open(filename, "w") as outfile:
            outfile.write("%g\n" % k)
            outfile.write("%g\n" % r_k)
            outfile.write("%g\n" % n_1)
            outfile.write("%g\n" % nu_k)
            outfile.write("%g\n" % integral_absolute)
            for i in range(n_1):
                values = self.point_array[-n_1+i]
                for value in values:
                    outfile.write("%g " % value)
                outfile.write("\n")
