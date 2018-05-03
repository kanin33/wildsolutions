import numpy as np
from Sequence import Sequence
from Dist_dU import Dist_dU
import time
import random

def generate_sequence(subsolution_function,
                      r_values,
                      t_min, t_max,
                      x_1_min=-1,
                      x_1_max=1,
                      max_iterations=100,
                      dist_filename="dist.txt",
                      dirname="test_dir",
                      z_k=None,
                      check_steps=True,
                      statistics_list=None):
    """Generate the sequence.

    Arguments:
    subsolution_function:
        return v_m, v_p, mu, delta, v_12, v_tilde, u_tilde, C_1, dist
    r_values:
        array of r_values to use in each iteration
    t_min, t_max:
        float
    x_1_min, x_1_max:
        float
    max_iterations:
        int, max iterations
    dist_filename:
        string, name on file where distance to boundary i stored
    dirname:
        string, name on directory to store solution
    z_k:
        Sequence object, if None, then a new is created
    check_steps:
        boolean, if False, do not check for any conditions in the construction
    statistics_list:
        list, used to store statistics on where the algorithm fails


    Returns:
    nu_k_values, list of nu-values corresponding to each iterations
    z_k, Sequence object representing the computed solution
    """
    subsolution_values = subsolution_function()
    mu_0     = subsolution_values[2][0]
    mu_1     = subsolution_values[2][1]
    v_tilde  = subsolution_values[5]
    u_tilde  = subsolution_values[6]
    C        = subsolution_values[7]

    # if not a Sequence object is provided, construct a new
    if z_k==None:
        dist_obj = Dist_dU(dist_filename) # dist_obj(x, y, z, q) = dist((x, y, z, q), dU)
        P        = [t_min, t_max, mu_0, mu_1, x_1_min, x_1_max]
        z_k      = Sequence(C, P, dist_obj, v_tilde, u_tilde, dirname=dirname)
    nu_k     = 0.5 # 2^-1
    nu_k_values = [nu_k]
    R_k      = 0

    accepted_iterations = 0


    # store statistics on what is failing
    # index 0 is inequality 1
    # index 1 is epsilon
    # index 2 is inequality 5
    # index 3 is success
    if statistics_list == None:
        statistics_list = [0, 0, 0, 0]
    for k in xrange(max_iterations):
        print "k: %d" % k
        # step 1 to 5
        r = r_values[k]
        print "r:", r
        accepted, R_k = advance_sequence(k,
                                         r,
                                         v_tilde,
                                         u_tilde,
                                         z_k,
                                         R_k,
                                         nu_k,
                                         statistics_list=statistics_list,
                                         check_steps=check_steps,
                                         offset_x1=random.random()*0.5,
                                         offset_x2=random.random()*0.5,
                                         prev_r=r_values[k-1] if k > 0 else 0)
        # step 6
        if accepted:
            print "computing next nu_k"
            nu_k = get_nu_k(accepted_iterations, z_k, n=10)
            print "nu = %f" % nu_k
            nu_k_values.append(nu_k)
            accepted_iterations += 1
            # write the accepted step to file
            z_k.write_step_to_file(nu_k)


    return nu_k_values, z_k

def advance_sequence(k,
                     r,
                     v_tilde,
                     u_tilde,
                     z_k,
                     R_k,
                     nu_k,
                     statistics_list,
                     check_steps=True,
                     offset_x1=0,
                     offset_x2=0,
                     prev_r=0,
                     ):
    """
    Advance the sequence one step. If it fails, return.


    Arguments:
    k:
        int, iteration number
    r:
        double, radius
    v_tilde:
        2-dim array
    u_tilde:
        2x2-dim array
    z_k:
        Sequence object
    R_k:
        double, current R_k
    nu_k:
        double, last nu_k
    check_steps:
        boolean, if False do not check anything
    offset_x1:
        double, how far from x1_min to start making cylinders
    offset_x2:
        double, how far from x2_min to start making cylinders
    prev_r:
        double, previous radius
    """
    # Timing
    start = time.time()

    # step 1
    print "Step 1"
    # alternate between two different ways to structure the points
    if k % 2 == 1 and len(z_k.r_values) > 0:
        prev_r = z_k.r_values[-1]
        points = find_points(z_k, prev_r, offset_x1=0, offset_x2=0)
        points = np.append(points,
                           find_points(z_k, prev_r,
                                       offset_x1=prev_r,
                                       offset_x2=prev_r),
                           axis=0)
        r = np.sqrt(2)/2*prev_r
    else:
        points = find_points(z_k, r, offset_x1=0, offset_x2=0)

    # if radius is too large
    if len(points) == 0:
        return False, R_k
    print "r: %.3f" % r

    # If not possible to find suitable points, change r, and redo step 1
    if check_steps and not inequality_step_1_b(r, z_k, points):
        print "Ineq 1 false"
        # GOTO step 1
        stop = time.time()
        print "Time used at step 1, iteration %d: %f" % (k, stop - start)
        statistics_list[0] += 1
        return False, R_k

    stop = time.time()
    print "Time used at step 1, iteration %d: %f" % (k, stop - start)

    start = time.time()
    # step 2 and 3
    print "Step 2 and 3"

    if check_steps:
        for point in points:
            d_j = dist_from_dU(z_k, point)
            epsilon_j = d_j - R_k
            if epsilon_j < 0:
                print "eps_j: %g" % epsilon_j
                print "d_j: %g" % d_j
                print "R_k: %g" % R_k
                # GOTO step 1
                stop = time.time()
                print "Time used at step 2 and 3, iteration %d: %f" % (k, stop - start)
                statistics_list[1] += 1
                return False, R_k
    print
    stop = time.time()
    print "Time used at step 2 and 3, iteration %d: %f" % (k, stop - start)

    # step 4
    start = time.time()
    print "Step 4"
    # Add the new points to the sequence
    z_k.add_points(points, r, n_trials=k, compute_integrals=check_steps)
    stop = time.time()
    print "Time used at step 4, iteration %d: %f" % (k, stop - start)

    start = time.time()
    # step 5
    print "Step 5"
    # if inequalities does not hold, back to step 1
    #TODO Check if the k value is correct
    if check_steps and not step_5_inequalities(z_k, nu_k):
        print "step5 false"
        # GOTO step 1
        stop = time.time()
        print "Time used at step 5, iteration %d: %f" % (k, stop - start)
        z_k.remove_last_points()
        statistics_list[2] += 1
        return False, R_k

    stop = time.time()
    print "Time used at step 5, iteration %d: %f" % (k, stop - start)

    # compute R_k+1
    print "Finding R_k"
    if check_steps:
        L_k = find_L_k(z_k, points, n=10)/20.
    else:
        L_k = 0
    print "L_k: %.4f, L_k*r: %.4f, r: %.4f" % (L_k, L_k*r, r)
    R_k = L_k*2*np.sqrt(2)*r
    stop = time.time()
    print "Time used to find L=%.2f, iteration %d: %f" % (L_k, k, stop - start)

    statistics_list[3] += 1
    return True, R_k

def find_points(z_k, r, offset_x1=0, offset_x2=0):
    """
    Compute finitely many points (x_j, t_j)

    The points are such that (t_j - r, t_j + r) x B_r(x_j) are disjoint
    and the inequality in step 1b holds.

    # The returned object is organized as a three dimentional array
    An empty list has boolean value False
    """
    # Assume nu_0 < 0, nu_1 > 0
    t_min   = z_k.t_min
    t_max   = z_k.t_max
    mu_0    = z_k.mu_0
    mu_1    = z_k.mu_1
    x_1_min = z_k.x_1_min
    x_1_max = z_k.x_1_max

    h = t_max - t_min         # height of set in t direction
    n_t = int(h/(2*r))        # number of cylinders in t direction

    # Note that the number of cylinders in x_1 direction
    # is constant in t, while the number of cylinders
    # in x_2 direction varies in t

    l_x_1 = x_1_max - x_1_min - offset_x1 # length of set in x_1 direction
    n_x_1 = int(l_x_1/(2*r))           # number of cylinders in x_1 direction

    # Start at t = t_min and compute the points for the first row
    # then compute the points for each row

    points = []
    for i in xrange(n_t):
        t_i     = t_min + i*2*r
        x_2_min = mu_0*t_i
        x_2_max = mu_1*t_i
        l_x_2   = x_2_max - x_2_min - offset_x1 if abs(x_2_max) > abs(x_2_min + offset_x1) else 0
        n_x_2   = int(l_x_2/(2*r))
        # x_1 direction
        for j in range(n_x_1):
            row_points = x_2_min + offset_x2 + r + np.arange(n_x_2)*2*r
            for point in row_points:
                # (t, x1, x2)
                points.append((t_i + r, x_1_min + offset_x1 + (2*j+1)*r, point))
    return np.array(points)

def get_nu_k(k, z_k, n=100):
    """Compute mu_k

    nu_k should satisfy:
        nu_k < 2^-k
        ||z_k - z_k*w_k|| < 2^-k
    """
    alpha_n = 1./3 # 4pi from C makes 4pi go awat
    const = np.sqrt(z_k.P_size)*alpha_n/0.035 # size of C
    t_min  = z_k.t_min;           t_max = z_k.t_max
    x1_min = z_k.x_1_min;        x1_max = z_k.x_1_max
    x2_min = z_k.t_max*z_k.mu_0; x2_max = z_k.t_max*z_k.mu_1
    t, X, Y = np.mgrid[t_min:t_max:n*1j, x1_min:x1_max:n*1j, x2_min:x2_max:n*1j]
    return 1./(2**k*const*z_k.max_grad(t, X, Y))

def find_L_k(z_k, points, n=100):
    # Compute L_k = max_grad
    t_min   = z_k.t_min
    t_max   = z_k.t_max
    x1_min  = z_k.x_1_min
    x1_max  = z_k.x_1_max
    mu_0    = z_k.mu_0
    mu_1    = z_k.mu_1
    x2_min  = t_max*mu_0
    x2_max  = t_max*mu_1
    t, X, Y = np.mgrid[t_min:t_max:1j*n,x1_min:x1_max:1j*n, x2_min:x2_max:1j*n]
    return z_k.max_grad(t, X, Y)

def step_5_inequalities(z_k, nu_k):
    """
    Check if the inequalitites in step 5 holds
    """

    c1 = 1./(400*np.sqrt(z_k.C)); c2 = 2./(4*np.pi) # todo: should not multiply by 2
    P_size = z_k.P_size

    lhs = z_k.integral_absolute()

    rhs = z_k.integral_absolute_prev()
    rhs += c1**2*c2**2/(2*P_size)*(z_k.C*P_size - z_k.integral_absolute_prev())**2

    inequality1 = lhs > rhs

    print "Inequality 5a:"
    print "|v_tilde-v^*|_L2^2                                     = %g" % lhs
    print "|v_tilde-v^k|_L2^2 + beta(C|P| - |v_tilde+v^k|_L2^2)^2 = %g" % rhs

    # convolution inequality
    k = len(z_k.r_values)
    z_star_z_k_int = z_k.integral_absolute_last_step()
    inequality2 = z_star_z_k_int < 2**(-k) # todo: maybe multiply k by 2?
    print "ineq5.2: %g < %g" % (z_star_z_k_int, 2**(-2*k))

    return inequality1 and inequality2

def step_1_integral(z_k):
    """
    Compute the integral in the inequality in step 1
    """
    # int_P |v_tilde + v_k|^2 dx
    print "Integrating step 1"
    integral = z_k.integral_absolute()
    P_size   = z_k.P_size
    print "integral 1: %f" % integral
    # 1/4pi int_P (C - |v_tilde + v_k|^2)dx
    return 1/(4*np.pi)*(P_size*z_k.C - integral)

def inequality_step_1_b(r, z_k, points):
    """
    Check if the inequality in step 1b holds.
    """
    lhs = 0
    for point in points:
        lhs = lhs + z_k.C - ((z_k.v_tilde[0] + z_k(*point)[0][0])**2 \
                       + (z_k.v_tilde[1] + z_k(*point)[0][1])**2)

    lhs = 2*r**3*lhs # todo: should not multiply by 2
    rhs = step_1_integral(z_k)
    print "Inequality 1b:"
    print "2*r^3 sum(C-|v_tilde+v_k(t_j,x_j)|^2)      = %g" % lhs
    print "c_2*iint(C-|v_tilde+v_k(t_j,x_j)|^2)dxdt = %g" % rhs

    return lhs > rhs

def dist_from_dU(z_k, point):
    return z_k.dist_from_dU(*point)
