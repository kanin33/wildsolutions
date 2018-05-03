import numpy as np
import time
import matlab.engine
import random
import sys

mlab = matlab.engine.start_matlab()

"""
This file contains very many functions, which purpose is to compute
different calculations and derivatives used in the construction of the
solutions.

The first six functions compute other things needed in the construction.
"""

def integrate_matlab(v_tilde, point_array):
    """
    Integrate iiint sumsum v_pj(2vtilde + v_k) dx1 dx2 dt

    :param point_array: Each column corresponds to one point
                        [r, t, x, y, p0, p1, p2, p3, N, eta1, eta2, eta3]
    :return:
    """
    print "Integrating"
    point_array = point_array.T
    v_tilde     = matlab.double(list(v_tilde))
    vars_k      = matlab.double(point_array.tolist())
    t1          = time.time()
    try:
        integral = mlab.cylinder_integral2(vars_k, v_tilde) # jekwl
    except:
        mlab.eval('exception = MException.last;', nargout=0)
        print mlab.eval('getReport(exception)')
        sys.exit()
    t2 = time.time()
    print "Time used to integrate: %.3f" % (t2 - t1)
    return integral

def points_on_K(C, offset=0):
    """
    Find five points on K evenly distributed
    """
    r   = np.sqrt(C)
    t   = 2*np.pi*np.arange(5)/5. + offset
    v1  = r * np.cos(t)
    v2  = r * np.sin(t)
    u11 = C/2. * np.cos(2*t)
    u12 = C/2. * np.sin(2*t)
    return v1, v2, u11, u12

def get_convex_combination(v_star, u_star, C):
    # Find five points on K
    v1_b, v2_b, u11_b, u12_b = points_on_K(C, offset=random.random()*2*np.pi)

    A = np.array([v1_b, v2_b, u11_b, u12_b, np.ones(5)])
    b = np.array([v_star[0], v_star[1], u_star[0], u_star[1], 1])

    alpha = np.linalg.solve(A, b)
    v     = np.array([v1_b, v2_b]).transpose()
    u     = np.array([u11_b, u12_b]).transpose()
    return alpha, v, u

def find_a_b(v_star, u_star, C):
    """Return a, b and alpha_j"""
    # test values (C = 1)
    alpha, v, u = get_convex_combination(v_star, u_star, C)
    i = np.argmax(alpha)
    b = v[i]
    # j = argmax(alpha*|v - b|)
    vb = v - b
    j  = np.argmax(np.multiply(alpha,
                               np.sqrt(vb[:,0]**2 + vb[:,1]**2)))
    a  = v[j]
    alpha_j = alpha[j]

    return a, b, alpha_j

def compute_eta_p(v, u, C):
    a, b, alpha_j  = find_a_b(v, u, C)
    a1, a2 = a
    b1, b2 = b
    adotb = a1*b1 + a2*b2

    a_abs = np.sqrt(a1**2 + a2**2)
    b_abs = np.sqrt(b1**2 + b2**2)

    t0t   = 2./3

    eta1  = -1./(a_abs*b_abs + adotb)**t0t * (a1 + b1)
    eta2  = -1./(a_abs*b_abs + adotb)**t0t * (a2 + b2)
    eta3  = (a_abs*b_abs + adotb)**(1./3)
    eta   = np.array([eta1, eta2, eta3])

    lambda_ = 0.5*alpha_j

    p = [(a1 - b1)*lambda_, (a2 - b2)*lambda_,
         (a1**2 - b1**2)*lambda_, (a1*a2 - b1*b2)*lambda_]

    return eta, lambda_, p

def compute_N(eta, k):
    # > 10/|eta|
    eta0, eta1, eta2 = eta
    eta_abs          = np.sqrt(eta0**2 + eta1**2 + eta2**2)
    N                = 10./eta_abs + 1
    return N+k+8 # todo: + noe?

def cutoff(x):
    # Function that is 1 on [-0.5, 0.5] and 0 at -1 and 1
    if isinstance(x, (float, int)):
        x = max(0.5, min(abs(x), 7./8))
    else:
        n = x.shape
        x = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x), 7./8*np.ones(n)))
    return -(8*x-7)**5*(286720*x**4 - 519680*x**3 + 358080*x**2 -110600*x+12901)/19683.

def cutoff_dx(x):
    if isinstance(x, (float, int)):
        x = max(0.5, min(abs(x), 7./8))
    else:
        n = x.shape
        x = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x), 7./8*np.ones(n)))
    return -(143360/2187.)*(16*x**2 - 22*x + 7)**4

def cutoff_dxx(x):
    if isinstance(x, (float, int)):
        x = max(0.5, min(abs(x), 7./8))
    else:
        n = x.shape
        x = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x), 7./8*np.ones(n)))
    return -(8*143360/2187.)*(16*x**2 - 22*x + 7)**3.*(16*x - 11)

def cutoff_dxxx(x):
    return -524*(16*x - 11)*(96*x - 66)*(16*x**2 - 22*x + 7)**2 - 8384*(16*x**2 - 22*x + 7)**3

def cutoff_dxxxx(x):
    return -(9175040/729.)*(2*x - 1)*(8*x - 7)*(16*x - 11)*(448*x**2 - 616*x + 205)

def phi_cut(t, x1, x2):
    result = cutoff(np.abs(t))*cutoff(np.sqrt(x1**2 + x2**2))
    return result

def phi_derivative_t(t, x1, x2):
    return cutoff_dx(np.abs(t))*cutoff(np.sqrt(x1**2 + x2**2))

def phi_derivative_x1(t, x1, x2):
    x_length = np.sqrt(x1**2 + x2**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff(np.abs(t))*x1/x_length*cutoff_dx(x_length)

def grad_cos_phi_abs(t, x1, x2, r, t_j, x1_j, x2_j, p0, p1, p2, p3, N, eta1, eta2, eta3, *args):
    tshape  = t.shape
    x1shape = x1.shape
    x2shape = x2.shape
    t  = t.reshape(tshape + (1,))
    x1 = x1.reshape(x1shape + (1,))
    x2 = x2.reshape(x2shape + (1,))
    t_new  = (t  - t_j)/r
    x1_new = (x1 - x1_j)/r
    x2_new = (x2 - x2_j)/r

    N_dot_eta = N*(x1_new*eta1 + x2_new*eta2 + t_new*eta3)
    cos_eta   = np.cos(N_dot_eta)
    sin_eta   = np.sin(N_dot_eta)
    cut_t     = cutoff(t_new)
    x_length  = np.sqrt(x1_new**2 + x2_new**2)
    cut_x     = cutoff(x_length)
    cut_derivative_t = cutoff_dx(t_new)
    cut_derivative_x = cutoff_dx(x_length)

    grad_cos_phi = np.array([cut_x*(-N*eta3*sin_eta*cut_t + cos_eta*cut_derivative_t),
                                  cut_t*(-N*eta1*sin_eta*cut_x + cos_eta*x1_new/x_length*cut_derivative_x),
                                  cut_t*(-N*eta2*sin_eta*cut_x + cos_eta*x2_new/x_length*cut_derivative_x)])

    v1_grad = p0*grad_cos_phi
    v2_grad = p1*grad_cos_phi

    u11_grad = p2*grad_cos_phi
    u12_grad = p3*grad_cos_phi

    return v1_grad, v2_grad, u11_grad, u12_grad

def v_p_j_hat(t, x1, x2, r, t_j, x1_j, x2_j, p0, p1, p2, p3, N, eta1, eta2, eta3, *args):
    tshape  = t.shape
    x1shape = x1.shape
    x2shape = x2.shape
    t  = t.reshape(tshape + (1,))
    x1 = x1.reshape(x1shape + (1,))
    x2 = x2.reshape(x2shape + (1,))
    t_new  = (t  - t_j)/r
    x1_new = (x1 - x1_j)/r
    x2_new = (x2 - x2_j)/r

    N_dot_eta = N*(x1_new*eta1 + x2_new*eta2 + t_new*eta3)
    cos_phi   = np.cos(N_dot_eta) \
                *phi_cut(t_new, x1_new, x2_new)

    v1 = r*p0*cos_phi
    v2 = r*p1*cos_phi

    u11 = r*p2*cos_phi
    u12 = r*p3*cos_phi

    return np.array([[v1, v2], [u11, u12]])

def v_p_j(t, x1, x2, r, t_j, x1_j, x2_j, p0, p1, p2, p3, N, eta1, eta2, eta3, a1, a2, b1, b2, lambda_, ampl=1):
    tshape  = t.shape
    x1shape = x1.shape
    x2shape = x2.shape
    t  = t.reshape(tshape + (1,))
    x1 = x1.reshape(x1shape + (1,))
    x2 = x2.reshape(x2shape + (1,))
    t_new  = (t  - t_j)/r
    x1_new = (x1 - x1_j)/r
    x2_new = (x2 - x2_j)/r

    #fact = N/(N+10)*0.5*(a1*b2 - a2*b1)
    #fact = 0.5*(a1*b2 - a2*b1)
    #np.random.seed(N)
    #coeff = random.random()
    #fact = coeff*0.5*(a1*b2 - a2*b1)
    #fact = ampl*r*0.5*(a1*b2 - a2*b1)
    #fact = r*0.5*(a1*b2 - a2*b1)
    fact = ampl*0.5*(a1*b2 - a2*b1)

    v1 = fact*lambda_*cutoff(t_new)*(
            dxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)*phi(     x1_new, x2_new)
        +   dxx_psi( N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dy_phi(  x1_new, x2_new)
        + 2*dxy_psi( N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dx_phi(  x1_new, x2_new)
        + 2*dx_psi(  N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dxy_phi( x1_new, x2_new)
        +   dy_psi(  N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dxx_phi( x1_new, x2_new)
        +   fpsi(    N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dxxy_phi(x1_new, x2_new)
        + 3*dyy_psi( N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dy_phi(  x1_new, x2_new)
        + 3*dy_psi(  N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dyy_phi( x1_new, x2_new)
        +   fpsi(    N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dyyy_phi(x1_new, x2_new)
        + dyyy_psi(  N, eta1, eta2, eta3, x1_new, x2_new, t_new)*phi(     x1_new, x2_new))

    v2 = -fact*lambda_*cutoff(t_new)*(
            dxxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)*phi(     x1_new, x2_new)
        + 3*dxx_psi( N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dx_phi(  x1_new, x2_new)
        + 3*dx_psi(  N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dxx_phi( x1_new, x2_new)
        +   fpsi(    N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dxxx_phi(x1_new, x2_new)
        +   dxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)*phi(     x1_new, x2_new)
        +   dyy_psi( N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dx_phi(  x1_new, x2_new)
        + 2*dxy_psi( N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dy_phi(  x1_new, x2_new)
        + 2*dy_psi(  N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dxy_phi( x1_new, x2_new)
        +   dx_psi(  N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dyy_phi( x1_new, x2_new)
        +   fpsi(    N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dxyy_phi(x1_new, x2_new))

    u11 = 2*fact*lambda_*(
          dxyt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)*phi(    x1_new, x2_new)*cutoff(   t_new)
        + dxy_psi( N, eta1, eta2, eta3, x1_new, x2_new, t_new)*phi(    x1_new, x2_new)*cutoff_dx(t_new)
        + dxt_psi( N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dy_phi( x1_new, x2_new)*cutoff(   t_new)
        + dyt_psi( N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dx_phi( x1_new, x2_new)*cutoff(   t_new)
        + dx_psi(  N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dy_phi( x1_new, x2_new)*cutoff_dx(t_new)
        + dy_psi(  N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dx_phi( x1_new, x2_new)*cutoff_dx(t_new)
        + dt_psi(  N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dxy_phi(x1_new, x2_new)*cutoff(   t_new)
        + fpsi(    N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dxy_phi(x1_new, x2_new)*cutoff_dx(t_new))

    u12 = fact*lambda_*(
          dyyt_psi( N, eta1, eta2, eta3, x1_new, x2_new, t_new)*phi(    x1_new, x2_new)*cutoff(   t_new)
        + dyy_psi(  N, eta1, eta2, eta3, x1_new, x2_new, t_new)*phi(    x1_new, x2_new)*cutoff_dx(t_new)
        + 2*dyt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dy_phi( x1_new, x2_new)*cutoff(   t_new)
        + 2*dy_psi( N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dy_phi( x1_new, x2_new)*cutoff_dx(t_new)
        + dt_psi(   N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dyy_phi(x1_new, x2_new)*cutoff(   t_new)
        + fpsi(     N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dyy_phi(x1_new, x2_new)*cutoff_dx(t_new)
        - dxxt_psi( N, eta1, eta2, eta3, x1_new, x2_new, t_new)*phi(    x1_new, x2_new)*cutoff(   t_new)
        - dxx_psi(  N, eta1, eta2, eta3, x1_new, x2_new, t_new)*phi(    x1_new, x2_new)*cutoff_dx(t_new)
        - 2*dxt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dx_phi( x1_new, x2_new)*cutoff(   t_new)
        - 2*dx_psi( N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dx_phi( x1_new, x2_new)*cutoff_dx(t_new)
        - dt_psi(   N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dxx_phi(x1_new, x2_new)*cutoff(   t_new)
        - fpsi(     N, eta1, eta2, eta3, x1_new, x2_new, t_new)*dxx_phi(x1_new, x2_new)*cutoff_dx(t_new))

    return np.array([[-v1, -v2], [u11, u12]])

def grad_v_p_j_new(t, x1, x2, r, t_j, x1_j, x2_j, p0, p1, p2, p3, N, eta1, eta2, eta3, a1, a2, b1, b2, lambda_):
    tshape  = t.shape
    x1shape = x1.shape
    x2shape = x2.shape
    t  = t.reshape(tshape + (1,))
    x1 = x1.reshape(x1shape + (1,))
    x2 = x2.reshape(x2shape + (1,))
    t_new  = (t  - t_j)/r
    x1_new = (x1 - x1_j)/r
    x2_new = (x2 - x2_j)/r

    fact = r*0.5*(a1*b2 - a2*b1)

    v1_grad = 1./r*fact*lambda_*np.array(
                   [v1_dx(N, eta1, eta2, eta3, x1_new, x2_new, t_new),
                    v1_dy(N, eta1, eta2, eta3, x1_new, x2_new, t_new),
                    v1_dt(N, eta1, eta2, eta3, x1_new, x2_new, t_new)])

    v2_grad = -1./r*fact*lambda_*np.array(
                    [v2_dx(N, eta1, eta2, eta3, x1_new, x2_new, t_new),
                     v2_dy(N, eta1, eta2, eta3, x1_new, x2_new, t_new),
                     v2_dt(N, eta1, eta2, eta3, x1_new, x2_new, t_new)])

    u11_grad = 2./r*fact*lambda_*np.array(
                    [u11_dx(N, eta1, eta2, eta3, x1_new, x2_new, t_new),
                     u11_dy(N, eta1, eta2, eta3, x1_new, x2_new, t_new),
                     u11_dt(N, eta1, eta2, eta3, x1_new, x2_new, t_new)])

    u12_grad = 1./r*fact*lambda_*np.array(
                    [u12_dx(N, eta1, eta2, eta3, x1_new, x2_new, t_new),
                     u12_dy(N, eta1, eta2, eta3, x1_new, x2_new, t_new),
                     u12_dt(N, eta1, eta2, eta3, x1_new, x2_new, t_new)])

    return v1_grad, v2_grad, u11_grad, u12_grad

def grad_v_p_j_b(t, x1, x2, r, t_j, x1_j, x2_j, p0, p1, p2, p3, N, eta1, eta2, eta3, a1, a2, b1, b2, lambda_, ampl=1):
    #np.random.seed(N)
    #coeff = np.random.random()
    #fact = r*0.5*(a1*b2 - a2*b1)/N # TODO: need to write why i divide by N
    #fact = ampl*r*0.5*(a1*b2 - a2*b1)/N # TODO: need to write why i divide by N
    fact = ampl*0.5*(a1*b2 - a2*b1)/N # TODO: need to write why i divide by N
    #fact = N/(N+10)*0.5*(a1*b2 - a2*b1)/N # TODO: need to write why i divide by N
    v1_grad_b = np.abs(1./r*fact*lambda_)**2*(v1_dx_b(N, eta1, eta2, eta3)**2
                                              + v1_dy_b(N, eta1, eta2, eta3)**2
                                              + v1_dt_b(N, eta1, eta2, eta3)**2)

    v2_grad_b = np.abs(1./r*fact*lambda_)**2*(v2_dx_b(N, eta1, eta2, eta3)**2
                                              + v2_dy_b(N, eta1, eta2, eta3)**2
                                              + v2_dt_b(N, eta1, eta2, eta3)**2)

    u11_grad_b = np.abs(2./r*fact*lambda_)**2*(u11_dx_b(N, eta1, eta2, eta3)**2
                                               + u11_dy_b(N, eta1, eta2, eta3)**2
                                               + u11_dt_b(N, eta1, eta2, eta3)**2)

    u12_grad_b = np.abs(1./r*fact*lambda_)**2*(u12_dx_b(N, eta1, eta2, eta3)**2
                                               + u12_dy_b(N, eta1, eta2, eta3)**2
                                               + u12_dt_b(N, eta1, eta2, eta3)**2)
    return v1_grad_b, v2_grad_b, u11_grad_b, u12_grad_b

def v1_dx_b(N, eta1, eta2, eta3 ):
    # compute bound of derivative
    value = cutoff_b()* \
            (phi_b()
             * dxxxy_psi_b(N, eta1, eta2, eta3 )
             + phi_b()
             * dxyyy_psi_b(N, eta1, eta2, eta3 )
             + fpsi_b(N, eta1, eta2, eta3 )
             * dxxxy_phi_b()
             + fpsi_b(N, eta1, eta2, eta3 )
             * dxyyy_phi_b()
             + 3 * dx_phi_b()
             * dxxy_psi_b(N, eta1, eta2, eta3 )
             + dx_phi_b()
             * dyyy_psi_b(N, eta1, eta2, eta3 )
             + dy_phi_b()
             * dxxx_psi_b(N, eta1, eta2, eta3 )
             + 3 * dy_phi_b()
             * dxyy_psi_b(N, eta1, eta2, eta3 )
             + 3 * dx_psi_b(N, eta1, eta2, eta3 )
             * dxxy_phi_b()
             + dx_psi_b(N, eta1, eta2, eta3 )
             * dyyy_phi_b()
             + dy_psi_b(N, eta1, eta2, eta3 )
             * dxxx_phi_b()
             + 3 * dy_psi_b(N, eta1, eta2, eta3 )
             * dxyy_phi_b()
             + 3 * dxx_phi_b()
             * dxy_psi_b(N, eta1, eta2, eta3 )
             + 3 * dxy_phi_b()
             * dxx_psi_b(N, eta1, eta2, eta3 )
             + 3 * dxy_phi_b()
             * dyy_psi_b(N, eta1, eta2, eta3 )
             + 3 * dyy_phi_b()
             * dxy_psi_b(N, eta1, eta2, eta3 ))
    return value

def v1_dy_b(N, eta1, eta2, eta3 ):
    value = cutoff_b()* \
            (phi_b()
             * dxxyy_psi_b(N, eta1, eta2, eta3 )
             + phi_b()
             * dyyyy_psi_b(N, eta1, eta2, eta3 )
             + fpsi_b(N, eta1, eta2, eta3 )
             * dxxyy_phi_b()
             + fpsi_b(N, eta1, eta2, eta3 )
             * dyyyy_phi_b()
             + 2*dx_phi_b()
             * dxyy_psi_b(N, eta1, eta2, eta3 )
             + 2*dy_phi_b()
             * dxxy_psi_b(N, eta1, eta2, eta3 )
             + 4*dy_phi_b()
             * dyyy_psi_b(N, eta1, eta2, eta3 )
             + 2*dx_psi_b(N, eta1, eta2, eta3 )
             * dxyy_phi_b()
             + 2*dy_psi_b(N, eta1, eta2, eta3 )
             * dxxy_phi_b()
             + 4*dy_psi_b(N, eta1, eta2, eta3 )
             * dyyy_phi_b()
             + dxx_phi_b()
             * dyy_psi_b(N, eta1, eta2, eta3 )
             + 4* dxy_phi_b()
             * dxy_psi_b(N, eta1, eta2, eta3 )
             + dyy_phi_b()
             * dxx_psi_b(N, eta1, eta2, eta3 )
             + 6 * dyy_phi_b()
             * dyy_psi_b(N, eta1, eta2, eta3 ))
    return value

def v1_dt_b(N, eta1, eta2, eta3 ):
    value =         phi_b()*cutoff_b() \
                    * dtxxy_psi_b(N, eta1, eta2, eta3 ) \
                    + phi_b()*cutoff_b() \
                      * dtyyy_psi_b(N, eta1, eta2, eta3 ) \
                    + fpsi_b(N, eta1, eta2, eta3 ) \
                      * dxxy_phi_b()*cutoff_dx_b() \
                    + fpsi_b(N, eta1, eta2, eta3 ) \
                      * dyyy_phi_b()*cutoff_dx_b() \
                    + phi_b()*cutoff_dx_b() \
                      * dxxy_psi_b(N, eta1, eta2, eta3 ) \
                    + phi_b()*cutoff_dx_b() \
                      * dyyy_psi_b(N, eta1, eta2, eta3 ) \
                    + 2*dx_phi_b()*cutoff_b() \
                      * dtxy_psi_b(N, eta1, eta2, eta3 ) \
                    + dy_phi_b()*cutoff_b() \
                      * dtxx_psi_b(N, eta1, eta2, eta3 ) \
                    + 3*dy_phi_b()*cutoff_b() \
                      * dyyt_psi_b(N, eta1, eta2, eta3 ) \
                    + dt_psi_b(N, eta1, eta2, eta3 ) \
                      * dxxy_phi_b()*cutoff_b() \
                    + dt_psi_b(N, eta1, eta2, eta3 ) \
                      * dyyy_phi_b()*cutoff_b() \
                    + 2*dx_psi_b(N, eta1, eta2, eta3 ) \
                      * dxy_phi_b()*cutoff_dx_b() \
                    + dy_psi_b(N, eta1, eta2, eta3 ) \
                      * dxx_phi_b()*cutoff_dx_b() \
                    + 3*dy_psi_b(N, eta1, eta2, eta3 ) \
                      * dyy_phi_b()*cutoff_dx_b() \
                    + 2*dx_phi_b()*cutoff_dx_b() \
                      * dxy_psi_b(N, eta1, eta2, eta3 ) \
                    + dy_phi_b()*cutoff_dx_b() \
                      * dxx_psi_b(N, eta1, eta2, eta3 ) \
                    + 3*dy_phi_b()*cutoff_dx_b() \
                      * dyy_psi_b(N, eta1, eta2, eta3 ) \
                    + dxx_phi_b()*cutoff_b() \
                      * dty_psi_b(N, eta1, eta2, eta3 ) \
                    + 2*dxy_phi_b()*cutoff_b() \
                      * dtx_psi_b(N, eta1, eta2, eta3 ) \
                    + 3*dyy_phi_b()*cutoff_b() \
                      * dty_psi_b(N, eta1, eta2, eta3 )
    return value

def v2_dx_b(N, eta1, eta2, eta3 ):
    value = cutoff_b()* \
            (phi_b() \
             * dxxxx_psi_b(N, eta1, eta2, eta3 )
             + phi_b() \
             * dxxyy_psi_b(N, eta1, eta2, eta3 )
             + fpsi_b(N, eta1, eta2, eta3 ) \
             * dxxxx_phi_b() \
             + fpsi_b(N, eta1, eta2, eta3 ) \
             * dxxyy_phi_b() \
             + 4 * dx_phi_b() \
             * dxxx_psi_b(N, eta1, eta2, eta3 ) \
             + 2 * dx_phi_b() \
             * dxyy_psi_b(N, eta1, eta2, eta3 ) \
             + 2 * dy_phi_b() \
             * dxxy_psi_b(N, eta1, eta2, eta3 ) \
             + 4 * dx_psi_b(N, eta1, eta2, eta3 ) \
             * dxxx_phi_b() \
             + 2 * dx_psi_b(N, eta1, eta2, eta3 ) \
             * dxyy_phi_b() \
             + 2 * dy_psi_b(N, eta1, eta2, eta3 ) \
             * dxxy_phi_b() \
             + 6 * dxx_phi_b() \
             * dxx_psi_b(N, eta1, eta2, eta3 ) \
             + dxx_phi_b() \
             * dyy_psi_b(N, eta1, eta2, eta3 ) \
             + 4 * dxy_phi_b() \
             * dxy_psi_b(N, eta1, eta2, eta3 ) \
             + dyy_phi_b() \
             * dxx_psi_b(N, eta1, eta2, eta3 ))

    return value

def v2_dy_b(N, eta1, eta2, eta3 ):
    value = cutoff_b()* \
            (phi_b()
             * dxxxy_psi_b(N, eta1, eta2, eta3 )
             + phi_b()
             * dxyyy_psi_b(N, eta1, eta2, eta3 )
             + fpsi_b(N, eta1, eta2, eta3 )
             * dxxxy_phi_b()
             + fpsi_b(N, eta1, eta2, eta3 )
             * dxyyy_phi_b()
             + 3 * dx_phi_b()
             * dxxy_psi_b(N, eta1, eta2, eta3 )
             + dx_phi_b()
             * dyyy_psi_b(N, eta1, eta2, eta3 )
             + dy_phi_b()
             * dxxx_psi_b(N, eta1, eta2, eta3 )
             + 3 * dy_phi_b()
             * dxyy_psi_b(N, eta1, eta2, eta3 )
             + 3 * dx_psi_b(N, eta1, eta2, eta3 )
             * dxxy_phi_b()
             + dx_psi_b(N, eta1, eta2, eta3 )
             * dyyy_phi_b()
             + dy_psi_b(N, eta1, eta2, eta3 )
             * dxxx_phi_b()
             + 3 * dy_psi_b(N, eta1, eta2, eta3 )
             * dxyy_phi_b()
             + 3 * dxx_phi_b()
             * dxy_psi_b(N, eta1, eta2, eta3 )
             + 3 * dxy_phi_b()
             * dxx_psi_b(N, eta1, eta2, eta3 )
             + 3 * dxy_phi_b()
             * dyy_psi_b(N, eta1, eta2, eta3 )
             + 3 * dyy_phi_b()
             * dxy_psi_b(N, eta1, eta2, eta3 ))

    return value

def v2_dt_b(N, eta1, eta2, eta3 ):
    value =         phi_b()*cutoff_b() \
                    * dtxxx_psi_b(N, eta1, eta2, eta3 ) \
                    + phi_b()*cutoff_b() \
                      * dtxyy_psi_b(N, eta1, eta2, eta3 ) \
                    + fpsi_b(N, eta1, eta2, eta3 ) \
                      * dxxx_phi_b()*cutoff_dx_b() \
                    + fpsi_b(N, eta1, eta2, eta3 ) \
                      * dxyy_phi_b()*cutoff_dx_b() \
                    + phi_b()*cutoff_dx_b() \
                      * dxxx_psi_b(N, eta1, eta2, eta3 ) \
                    + phi_b()*cutoff_dx_b() \
                      * dxyy_psi_b(N, eta1, eta2, eta3 ) \
                    + 3*dx_phi_b()*cutoff_b() \
                      * dtxx_psi_b(N, eta1, eta2, eta3 ) \
                    + dx_phi_b()*cutoff_b() \
                      * dtyy_psi_b(N, eta1, eta2, eta3 ) \
                    + 2*dy_phi_b()*cutoff_b() \
                      * dtxy_psi_b(N, eta1, eta2, eta3 ) \
                    + dt_psi_b(N, eta1, eta2, eta3 ) \
                      * dxxx_phi_b()*cutoff_b() \
                    + dt_psi_b(N, eta1, eta2, eta3 ) \
                      * dxyy_phi_b()*cutoff_b() \
                    + 3*dx_psi_b(N, eta1, eta2, eta3 ) \
                      * dxx_phi_b()*cutoff_dx_b() \
                    + dx_psi_b(N, eta1, eta2, eta3 ) \
                      * dyy_phi_b()*cutoff_dx_b() \
                    + 2*dy_psi_b(N, eta1, eta2, eta3 ) \
                      * dxy_phi_b()*cutoff_dx_b() \
                    + 3*dx_phi_b()*cutoff_dx_b() \
                      * dxx_psi_b(N, eta1, eta2, eta3 ) \
                    + dx_phi_b()*cutoff_dx_b() \
                      * dyy_psi_b(N, eta1, eta2, eta3 ) \
                    + 2*dy_phi_b()*cutoff_dx_b() \
                      * dxy_psi_b(N, eta1, eta2, eta3 ) \
                    + 3*dxx_phi_b()*cutoff_b() \
                      * dtx_psi_b(N, eta1, eta2, eta3 ) \
                    + 2*dxy_phi_b()*cutoff_b() \
                      * dty_psi_b(N, eta1, eta2, eta3 ) \
                    + dyy_phi_b()*cutoff_b() \
                      * dtx_psi_b(N, eta1, eta2, eta3 )

    return value

def u11_dx_b(N, eta1, eta2, eta3 ):
    value =         phi_b()*cutoff_b() \
                    * dtxxy_psi_b(N, eta1, eta2, eta3 ) \
                    + fpsi_b(N, eta1, eta2, eta3 ) \
                      * dxxy_phi_b()*cutoff_dx_b() \
                    + phi_b()*cutoff_dx_b() \
                      * dxxy_psi_b(N, eta1, eta2, eta3 ) \
                    + 2*dx_phi_b()*cutoff_b() \
                      * dtxy_psi_b(N, eta1, eta2, eta3 ) \
                    + dy_phi_b()*cutoff_b() \
                      * dtxx_psi_b(N, eta1, eta2, eta3 ) \
                    + dt_psi_b(N, eta1, eta2, eta3 ) \
                      * dxxy_phi_b()*cutoff_b() \
                    + 2*dx_psi_b(N, eta1, eta2, eta3 ) \
                      * dxy_phi_b()*cutoff_dx_b() \
                    + dy_psi_b(N, eta1, eta2, eta3 ) \
                      * dxx_phi_b()*cutoff_dx_b() \
                    + 2*dx_phi_b()*cutoff_dx_b() \
                      * dxy_psi_b(N, eta1, eta2, eta3 ) \
                    + dy_phi_b()*cutoff_dx_b() \
                      * dxx_psi_b(N, eta1, eta2, eta3 ) \
                    + dxx_phi_b()*cutoff_b() \
                      * dty_psi_b(N, eta1, eta2, eta3 ) \
                    + 2*dxy_phi_b()*cutoff_b() \
                      * dtx_psi_b(N, eta1, eta2, eta3 )

    return value

def u11_dy_b(N, eta1, eta2, eta3 ):
    value =         phi_b()*cutoff_b() \
                    * dtxyy_psi_b(N, eta1, eta2, eta3 ) \
                    + fpsi_b(N, eta1, eta2, eta3 ) \
                      * dxyy_phi_b()*cutoff_dx_b() \
                    + phi_b()*cutoff_dx_b() \
                      * dxyy_psi_b(N, eta1, eta2, eta3 ) \
                    + dx_phi_b()*cutoff_b() \
                      * dtyy_psi_b(N, eta1, eta2, eta3 ) \
                    + 2 * dy_phi_b()*cutoff_b() \
                      * dtxy_psi_b(N, eta1, eta2, eta3 ) \
                    + dt_psi_b(N, eta1, eta2, eta3 ) \
                      * dxyy_phi_b()*cutoff_b() \
                    + dx_psi_b(N, eta1, eta2, eta3 ) \
                      * dyy_phi_b()*cutoff_dx_b() \
                    + 2 * dy_psi_b(N, eta1, eta2, eta3 ) \
                      * dxy_phi_b()*cutoff_dx_b() \
                    + dx_phi_b()*cutoff_dx_b() \
                      * dyy_psi_b(N, eta1, eta2, eta3 ) \
                    + 2 * dy_phi_b()*cutoff_dx_b() \
                      * dxy_psi_b(N, eta1, eta2, eta3 ) \
                    + 2 * dxy_phi_b()*cutoff_b() \
                      * dty_psi_b(N, eta1, eta2, eta3 ) \
                    + dyy_phi_b()*cutoff_b() \
                      * dtx_psi_b(N, eta1, eta2, eta3 )

    return value

def u11_dt_b(N, eta1, eta2, eta3 ):
    value =         phi_b()*cutoff_b() \
                    * dttxy_psi_b(N, eta1, eta2, eta3 ) \
                    + fpsi_b(N, eta1, eta2, eta3 ) \
                      * dxy_phi_b()*cutoff_dxx_b() \
                    + 2*phi_b()*cutoff_dx_b() \
                      * dtxy_psi_b(N, eta1, eta2, eta3 ) \
                    + dx_phi_b()*cutoff_b() \
                      * dtty_psi_b(N, eta1, eta2, eta3 ) \
                    + dy_phi_b()*cutoff_b() \
                      * dttx_psi_b(N, eta1, eta2, eta3 ) \
                    + 2*dt_psi_b(N, eta1, eta2, eta3 ) \
                      * dxy_phi_b()*cutoff_dx_b() \
                    + dx_psi_b(N, eta1, eta2, eta3 ) \
                      * dy_phi_b()*cutoff_dxx_b() \
                    + dy_psi_b(N, eta1, eta2, eta3 ) \
                      * dx_phi_b()*cutoff_dxx_b() \
                    + phi_b()*cutoff_dxx_b() \
                      * dxy_psi_b(N, eta1, eta2, eta3 ) \
                    + 2*dx_phi_b()*cutoff_dx_b() \
                      * dty_psi_b(N, eta1, eta2, eta3 ) \
                    + 2*dy_phi_b()*cutoff_dx_b() \
                      * dtx_psi_b(N, eta1, eta2, eta3 ) \
                    + dxy_phi_b()*cutoff_b() \
                      * dtt_psi_b(N, eta1, eta2, eta3 )

    return value

def u12_dx_b(N, eta1, eta2, eta3 ):
    value =        -phi_b()*cutoff_b() \
                   * dtxxx_psi_b(N, eta1, eta2, eta3 ) \
                   + phi_b()*cutoff_b() \
                     * dtxyy_psi_b(N, eta1, eta2, eta3 ) \
                   - fpsi_b(N, eta1, eta2, eta3 ) \
                     * dxxx_phi_b()*cutoff_dx_b() \
                   + fpsi_b(N, eta1, eta2, eta3 ) \
                     * dxyy_phi_b()*cutoff_dx_b() \
                   - phi_b()*cutoff_dx_b() \
                     * dxxx_psi_b(N, eta1, eta2, eta3 ) \
                   + phi_b()*cutoff_dx_b() \
                     * dxyy_psi_b(N, eta1, eta2, eta3 ) \
                   - 3*dx_phi_b()*cutoff_b() \
                     * dtxx_psi_b(N, eta1, eta2, eta3 ) \
                   + dx_phi_b()*cutoff_b() \
                     * dtyy_psi_b(N, eta1, eta2, eta3 ) \
                   + 2*dy_phi_b()*cutoff_b() \
                     * dtxy_psi_b(N, eta1, eta2, eta3 ) \
                   - dt_psi_b(N, eta1, eta2, eta3 ) \
                     * dxxx_phi_b()*cutoff_b() \
                   + dt_psi_b(N, eta1, eta2, eta3 ) \
                     * dxyy_phi_b()*cutoff_b() \
                   - 3*dx_psi_b(N, eta1, eta2, eta3 ) \
                     * dxx_phi_b()*cutoff_dx_b() \
                   + dx_psi_b(N, eta1, eta2, eta3 ) \
                     * dyy_phi_b()*cutoff_dx_b() \
                   + 2*dy_psi_b(N, eta1, eta2, eta3 ) \
                     * dxy_phi_b()*cutoff_dx_b() \
                   - 3*dx_phi_b()*cutoff_dx_b() \
                     * dxx_psi_b(N, eta1, eta2, eta3 ) \
                   + dx_phi_b()*cutoff_dx_b() \
                     * dyy_psi_b(N, eta1, eta2, eta3 ) \
                   + 2*dy_phi_b()*cutoff_dx_b() \
                     * dxy_psi_b(N, eta1, eta2, eta3 ) \
                   - 3*dxx_phi_b()*cutoff_b() \
                     * dtx_psi_b(N, eta1, eta2, eta3 ) \
                   + 2*dxy_phi_b()*cutoff_b() \
                     * dty_psi_b(N, eta1, eta2, eta3 ) \
                   + dyy_phi_b()*cutoff_b() \
                     *dtx_psi_b(N, eta1, eta2, eta3 )

    return value

def u12_dy_b(N, eta1, eta2, eta3 ):
    value =        -phi_b()*cutoff_b() \
                   * dtxxy_psi_b(N, eta1, eta2, eta3 ) \
                   + phi_b()*cutoff_b() \
                     * dtyyy_psi_b(N, eta1, eta2, eta3 ) \
                   - fpsi_b(N, eta1, eta2, eta3 ) \
                     * dxxy_phi_b()*cutoff_dx_b() \
                   + fpsi_b(N, eta1, eta2, eta3 ) \
                     * dyyy_phi_b()*cutoff_dx_b() \
                   - phi_b()*cutoff_dx_b() \
                     * dxxy_psi_b(N, eta1, eta2, eta3 ) \
                   + phi_b()*cutoff_dx_b() \
                     * dyyy_psi_b(N, eta1, eta2, eta3 ) \
                   - 2*dx_phi_b()*cutoff_b() \
                     * dtxy_psi_b(N, eta1, eta2, eta3 ) \
                   - dy_phi_b()*cutoff_b() \
                     * dtxx_psi_b(N, eta1, eta2, eta3 ) \
                   + 3*dy_phi_b()*cutoff_b() \
                     * dtyy_psi_b(N, eta1, eta2, eta3 ) \
                   - dt_psi_b(N, eta1, eta2, eta3 ) \
                     * dxxy_phi_b()*cutoff_b() \
                   + dt_psi_b(N, eta1, eta2, eta3 ) \
                     * dyyy_phi_b()*cutoff_b() \
                   - 2*dx_psi_b(N, eta1, eta2, eta3 ) \
                     * dxy_phi_b()*cutoff_dx_b() \
                   - dy_psi_b(N, eta1, eta2, eta3 ) \
                     * dxx_phi_b()*cutoff_dx_b() \
                   + 3*dy_psi_b(N, eta1, eta2, eta3 ) \
                     * dyy_phi_b()*cutoff_dx_b() \
                   - 2*dx_phi_b()*cutoff_dx_b() \
                     * dxy_psi_b(N, eta1, eta2, eta3 ) \
                   - dy_phi_b()*cutoff_dx_b() \
                     * dxx_psi_b(N, eta1, eta2, eta3 ) \
                   + 3*dy_phi_b()*cutoff_dx_b() \
                     * dyy_psi_b(N, eta1, eta2, eta3 ) \
                   - dxx_phi_b()*cutoff_b() \
                     * dty_psi_b(N, eta1, eta2, eta3 ) \
                   - 2*dxy_phi_b()*cutoff_b() \
                     * dtx_psi_b(N, eta1, eta2, eta3 ) \
                   + 3*dyy_phi_b()*cutoff_b() \
                     * dty_psi_b(N, eta1, eta2, eta3 )

    return value

def u12_dt_b(N, eta1, eta2, eta3 ):
    value =        -phi_b()*cutoff_b() \
                   * dttxx_psi_b(N, eta1, eta2, eta3 ) \
                   + phi_b()*cutoff_b() \
                     * dttyy_psi_b(N, eta1, eta2, eta3 ) \
                   - fpsi_b(N, eta1, eta2, eta3 ) \
                     * dxx_phi_b()*cutoff_dxx_b() \
                   + fpsi_b(N, eta1, eta2, eta3 ) \
                     * dyy_phi_b()*cutoff_dxx_b() \
                   - 2*phi_b()*cutoff_dx_b() \
                     * dtxx_psi_b(N, eta1, eta2, eta3 ) \
                   + 2*phi_b()*cutoff_dx_b() \
                     * dtyy_psi_b(N, eta1, eta2, eta3 ) \
                   - 2*dx_phi_b()*cutoff_b() \
                     * dttx_psi_b(N, eta1, eta2, eta3 ) \
                   + 2*dy_phi_b()*cutoff_b() \
                     * dtty_psi_b(N, eta1, eta2, eta3 ) \
                   - 2*dt_psi_b(N, eta1, eta2, eta3 ) \
                     * dxx_phi_b()*cutoff_dx_b() \
                   + 2*dt_psi_b(N, eta1, eta2, eta3 ) \
                     * dyy_phi_b()*cutoff_dx_b() \
                   - 2*dx_psi_b(N, eta1, eta2, eta3 ) \
                     * dx_phi_b()*cutoff_dxx_b() \
                   + 2*dy_psi_b(N, eta1, eta2, eta3 ) \
                     * dy_phi_b()*cutoff_dxx_b() \
                   - phi_b()*cutoff_dxx_b() \
                     * dxx_psi_b(N, eta1, eta2, eta3 ) \
                   + phi_b()*cutoff_dxx_b() \
                     * dyy_psi_b(N, eta1, eta2, eta3 ) \
                   - 4*dx_phi_b()*cutoff_dx_b() \
                     * dtx_psi_b(N, eta1, eta2, eta3 ) \
                   + 4*dy_phi_b()*cutoff_dx_b() \
                     * dty_psi_b(N, eta1, eta2, eta3 ) \
                   - dxx_phi_b()*cutoff_b() \
                     * dtt_psi_b(N, eta1, eta2, eta3 ) \
                   + dyy_phi_b()*cutoff_b() \
                     * dtt_psi_b(N, eta1, eta2, eta3 )

    return value

def v1_dx(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    value = cutoff(t_new)*\
                   (phi(x1_new, x2_new)
            * dxxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
                  + phi(x1_new, x2_new)
            * dxyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
            * dxxxy_phi(x1_new, x2_new)
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
            * dxyyy_phi(x1_new, x2_new)
           + 3 * dx_phi(x1_new, x2_new)
             * dxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
               + dx_phi(x1_new, x2_new)
             * dyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
               + dy_phi(x1_new, x2_new)
             * dxxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
           + 3 * dy_phi(x1_new, x2_new)
             * dxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
           + 3 * dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             * dxxy_phi(x1_new, x2_new)
               + dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             * dyyy_phi(x1_new, x2_new)
               + dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             * dxxx_phi(x1_new, x2_new)
           + 3 * dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             * dxyy_phi(x1_new, x2_new)
          + 3 * dxx_phi(x1_new, x2_new)
              * dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
          + 3 * dxy_phi(x1_new, x2_new)
              * dxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
          + 3 * dxy_phi(x1_new, x2_new)
              * dyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
          + 3 * dyy_phi(x1_new, x2_new)
              * dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new))
    return value

def v1_dy(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    value = cutoff(t_new)* \
                   (phi(x1_new, x2_new)
            * dxxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
                  + phi(x1_new, x2_new)
            * dyyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
            * dxxyy_phi(x1_new, x2_new)
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
            * dyyyy_phi(x1_new, x2_new)
             + 2*dx_phi(x1_new, x2_new)
             * dxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             + 2*dy_phi(x1_new, x2_new)
             * dxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             + 4*dy_phi(x1_new, x2_new)
             * dyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             + 2*dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             * dxyy_phi(x1_new, x2_new)
             + 2*dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             * dxxy_phi(x1_new, x2_new)
             + 4*dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             * dyyy_phi(x1_new, x2_new)
              + dxx_phi(x1_new, x2_new)
              * dyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
           + 4* dxy_phi(x1_new, x2_new)
              * dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
              + dyy_phi(x1_new, x2_new)
              * dxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
          + 6 * dyy_phi(x1_new, x2_new)
              * dyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new))
    return value

def v1_dt(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    value =         phi(x1_new, x2_new)*cutoff(t_new)\
            * dtxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                  + phi(x1_new, x2_new)*cutoff(t_new)\
            * dtyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxxy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dyyy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
                  + phi(x1_new, x2_new)*cutoff_dx(t_new)\
             * dxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                  + phi(x1_new, x2_new)*cutoff_dx(t_new)\
             * dyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 2*dx_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dy_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 3*dy_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxxy_phi(x1_new, x2_new)*cutoff(t_new)\
               + dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dyyy_phi(x1_new, x2_new)*cutoff(t_new)\
             + 2*dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
               + dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
             + 3*dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dyy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
             + 2*dx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 3*dy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              + dxx_phi(x1_new, x2_new)*cutoff(t_new)\
              * dty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
            + 2*dxy_phi(x1_new, x2_new)*cutoff(t_new)\
              * dtx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
            + 3*dyy_phi(x1_new, x2_new)*cutoff(t_new)\
              * dty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
    return value

def v2_dx(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    value = cutoff(t_new)* \
                   (phi(x1_new, x2_new)\
            * dxxxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
                 + phi(x1_new, x2_new)\
           * dxxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
           * dxxxx_phi(x1_new, x2_new)\
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
           * dxxyy_phi(x1_new, x2_new)\
          + 4 * dx_phi(x1_new, x2_new)\
            * dxxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
          + 2 * dx_phi(x1_new, x2_new)\
            * dxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
          + 2 * dy_phi(x1_new, x2_new)\
            * dxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
          + 4 * dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
            * dxxx_phi(x1_new, x2_new)\
          + 2 * dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
            * dxyy_phi(x1_new, x2_new)\
          + 2 * dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
            * dxxy_phi(x1_new, x2_new)\
         + 6 * dxx_phi(x1_new, x2_new)\
             * dxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + dxx_phi(x1_new, x2_new)\
             * dyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
         + 4 * dxy_phi(x1_new, x2_new)\
             * dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + dyy_phi(x1_new, x2_new)\
             * dxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new))

    return value

def v2_dy(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    value = cutoff(t_new)* \
                   (phi(x1_new, x2_new)
            * dxxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
                  + phi(x1_new, x2_new)
            * dxyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
            * dxxxy_phi(x1_new, x2_new)
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
            * dxyyy_phi(x1_new, x2_new)
           + 3 * dx_phi(x1_new, x2_new)
             * dxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
               + dx_phi(x1_new, x2_new)
             * dyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
               + dy_phi(x1_new, x2_new)
             * dxxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
           + 3 * dy_phi(x1_new, x2_new)
             * dxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
           + 3 * dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             * dxxy_phi(x1_new, x2_new)
               + dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             * dyyy_phi(x1_new, x2_new)
               + dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             * dxxx_phi(x1_new, x2_new)
           + 3 * dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
             * dxyy_phi(x1_new, x2_new)
          + 3 * dxx_phi(x1_new, x2_new)
              * dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
          + 3 * dxy_phi(x1_new, x2_new)
              * dxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
          + 3 * dxy_phi(x1_new, x2_new)
              * dyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
          + 3 * dyy_phi(x1_new, x2_new)
              * dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new))

    return value

def v2_dt(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    value =         phi(x1_new, x2_new)*cutoff(t_new)\
            * dtxxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                  + phi(x1_new, x2_new)*cutoff(t_new)\
            * dtxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxxx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxyy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
                  + phi(x1_new, x2_new)*cutoff_dx(t_new)\
             * dxxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                  + phi(x1_new, x2_new)*cutoff_dx(t_new)\
             * dxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 3*dx_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dx_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 2*dy_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxxx_phi(x1_new, x2_new)*cutoff(t_new)\
               + dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxyy_phi(x1_new, x2_new)*cutoff(t_new)\
             + 3*dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
               + dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dyy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
             + 2*dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
             + 3*dx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 2*dy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
            + 3*dxx_phi(x1_new, x2_new)*cutoff(t_new)\
              * dtx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
            + 2*dxy_phi(x1_new, x2_new)*cutoff(t_new)\
              * dty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              + dyy_phi(x1_new, x2_new)*cutoff(t_new)\
              * dtx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)

    return value

def u11_dx(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    value =         phi(x1_new, x2_new)*cutoff(t_new) \
            * dtxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new) \
             * dxxy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
                  + phi(x1_new, x2_new)*cutoff_dx(t_new)\
             * dxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 2*dx_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dy_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxxy_phi(x1_new, x2_new)*cutoff(t_new)\
             + 2*dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
               + dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
             + 2*dx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              + dxx_phi(x1_new, x2_new)*cutoff(t_new)\
              * dty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
            + 2*dxy_phi(x1_new, x2_new)*cutoff(t_new)\
              * dtx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)

    return value

def u11_dy(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    value =         phi(x1_new, x2_new)*cutoff(t_new)\
            * dtxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                  + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxyy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
                  + phi(x1_new, x2_new)*cutoff_dx(t_new)\
             * dxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dx_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
           + 2 * dy_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxyy_phi(x1_new, x2_new)*cutoff(t_new)\
               + dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dyy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
           + 2 * dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
               + dx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
           + 2 * dy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
          + 2 * dxy_phi(x1_new, x2_new)*cutoff(t_new)\
              * dty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              + dyy_phi(x1_new, x2_new)*cutoff(t_new)\
              * dtx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)

    return value

def u11_dt(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    value =         phi(x1_new, x2_new)*cutoff(t_new)\
            * dttxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new) \
                  + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxy_phi(x1_new, x2_new)*cutoff_dxx(t_new) \
                + 2*phi(x1_new, x2_new)*cutoff_dx(t_new)\
             * dtxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dx_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dy_phi(x1_new, x2_new)*cutoff(t_new)\
             * dttx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 2*dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
               + dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               * dy_phi(x1_new, x2_new)*cutoff_dxx(t_new)\
               + dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               * dx_phi(x1_new, x2_new)*cutoff_dxx(t_new)\
                  + phi(x1_new, x2_new)*cutoff_dxx(t_new)\
              * dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 2*dx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 2*dy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dtx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              + dxy_phi(x1_new, x2_new)*cutoff(t_new)\
              * dtt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)

    return value

def u12_dx(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    value =        -phi(x1_new, x2_new)*cutoff(t_new)\
            * dtxxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                  + phi(x1_new, x2_new)*cutoff(t_new)\
            * dtxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                 - fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxxx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxyy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
                  - phi(x1_new, x2_new)*cutoff_dx(t_new)\
             * dxxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                  + phi(x1_new, x2_new)*cutoff_dx(t_new)\
             * dxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             - 3*dx_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dx_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 2*dy_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               - dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxxx_phi(x1_new, x2_new)*cutoff(t_new)\
               + dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxyy_phi(x1_new, x2_new)*cutoff(t_new)\
             - 3*dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
               + dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dyy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
             + 2*dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
             - 3*dx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               + dx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 2*dy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
            - 3*dxx_phi(x1_new, x2_new)*cutoff(t_new)\
              * dtx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
            + 2*dxy_phi(x1_new, x2_new)*cutoff(t_new)\
              * dty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              + dyy_phi(x1_new, x2_new)*cutoff(t_new)\
               *dtx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)

    return value

def u12_dy(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    value =        -phi(x1_new, x2_new)*cutoff(t_new)\
            * dtxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                  + phi(x1_new, x2_new)*cutoff(t_new)\
            * dtyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                  - fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxxy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
                  + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dyyy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
                  - phi(x1_new, x2_new)*cutoff_dx(t_new)\
             * dxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                  + phi(x1_new, x2_new)*cutoff_dx(t_new)\
             * dyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             - 2*dx_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               - dy_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 3*dy_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               - dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dxxy_phi(x1_new, x2_new)*cutoff(t_new)\
               + dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             * dyyy_phi(x1_new, x2_new)*cutoff(t_new)\
             - 2*dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
               - dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
             + 3*dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dyy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
             - 2*dx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               - dy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 3*dy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              - dxx_phi(x1_new, x2_new)*cutoff(t_new)\
              * dty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
            - 2*dxy_phi(x1_new, x2_new)*cutoff(t_new)\
              * dtx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
            + 3*dyy_phi(x1_new, x2_new)*cutoff(t_new)\
              * dty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)

    return value

def u12_dt(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    value =        -phi(x1_new, x2_new)*cutoff(t_new)\
            * dttxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                  + phi(x1_new, x2_new)*cutoff(t_new)\
            * dttyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                 - fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxx_phi(x1_new, x2_new)*cutoff_dxx(t_new)\
                 + fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dyy_phi(x1_new, x2_new)*cutoff_dxx(t_new)\
                - 2*phi(x1_new, x2_new)*cutoff_dx(t_new)\
             * dtxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                + 2*phi(x1_new, x2_new)*cutoff_dx(t_new)\
             * dtyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             - 2*dx_phi(x1_new, x2_new)*cutoff(t_new)\
             * dttx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 2*dy_phi(x1_new, x2_new)*cutoff(t_new)\
             * dtty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             - 2*dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dxx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
             + 2*dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              * dyy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
             - 2*dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               * dx_phi(x1_new, x2_new)*cutoff_dxx(t_new)\
             + 2*dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
               * dy_phi(x1_new, x2_new)*cutoff_dxx(t_new)\
                  - phi(x1_new, x2_new)*cutoff_dxx(t_new)\
              * dxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
                  + phi(x1_new, x2_new)*cutoff_dxx(t_new)\
              * dyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             - 4*dx_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dtx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
             + 4*dy_phi(x1_new, x2_new)*cutoff_dx(t_new)\
              * dty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              - dxx_phi(x1_new, x2_new)*cutoff(t_new)\
              * dtt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)\
              + dyy_phi(x1_new, x2_new)*cutoff(t_new) \
              * dtt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)

    return value


def dyyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta2**4*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))*N

def dxxxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta1**4*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))*N

def dxxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta1**2*eta2**2*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))*N

def dxxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta1**3*eta2*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))*N

def dxyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta1*eta2**3*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))*N

def dtxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta3*eta1**2*eta2*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))*N

def dtxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta3*eta1*eta2**2*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))*N

def dtyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta3*eta2**3*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))*N

def dtxxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta3*eta1**3*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))*N

def dttxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta3**2*eta1*eta2*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))*N

def dttxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta3**2*eta1**2*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))*N

def dttyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta3**2*eta2**2*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))*N

def dtxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta3*eta1**2*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))

def dtyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta3*eta2**2*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))

def dtxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta3*eta1*eta2*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))

def dxxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta1**3*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))

def dyyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta2**3*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))

def dxxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta1**2*eta2*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))

def dxyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta1*eta2**2*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))

def dttx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta3**2*eta1*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))

def dtty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta3**2*eta2*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))

def dxx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta1**2*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))/N

def dyy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta2**2*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))/N

def dxy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta1*eta2*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))/N

def dtx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta3*eta1*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))/N

def dty_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta3*eta2*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))/N

def dx_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta1*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))/N**2

def dy_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta2*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))/N**2

def fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))/N**3

def dyyt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta2**2*eta3*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))

def dxyt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta1*eta2*eta3*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))

def dxxt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta1**2*eta3*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))

def dxt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta1*eta3*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))/N

def dyt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta2*eta3*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))/N

def dt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return -eta3*np.cos(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))/N**2

def dtt_psi(N, eta1, eta2, eta3, x1_new, x2_new, t_new):
    return eta3**2*np.sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))/N

def dxxyy_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dxxxx(x_length)*x**2*y**2/x_length**4 \
           + cutoff_dxxx(x_length)*(-4*x**2*y**2 + x**4 + y**4)/x_length**5 \
           + cutoff_dxx(x_length)*(7*x**2*y**2 - 4*x**4 - 4*y**2)/x_length**6 \
           + cutoff_dx(x_length)*(7*x**2*y**2 - 4*x**4 - 4*y**2)/x_length**7

def dxxxy_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dxxxx(x_length)*y*x**3/x_length**4 \
           + cutoff_dxxx(x_length)*y*x**3*(3*x**2 + y**2)/x_length**5 \
           + cutoff_dxx(x_length)*3*y*x*(2*x**2 - 3*y**2)/x_length**6 \
           + cutoff_dx(x_length)*y*x*3*(-2*x**2 + 2*y**2)/x_length**7

def dxyyy_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dxxxx(x_length)*x*y**3/x_length**4 \
           + cutoff_dxxx(x_length)*x*y**3*(3*y**2 + x**2)/x_length**5 \
           + cutoff_dxx(x_length)*3*x*y*(2*y**2 - 3*x**2)/x_length**6 \
           + cutoff_dx(x_length)*x*y*3*(-2*y**2 + 2*x**2)/x_length**7

def dyyyy_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dxxxx(x_length)*y**4/x_length**4 \
           + cutoff_dxxx(x_length)*6*y**2*x**2/x_length**5 \
           + cutoff_dxx(x_length)*(-18*y**2*x**2 + 3*x**4)/x_length**6 \
           + cutoff_dx(x_length)*(3*y**4 + 18*y**2*x**2)/x_length**7

def dxxxx_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dxxxx(x_length)*x**4/x_length**4 \
           + cutoff_dxxx(x_length)*6*x**2*y**2/x_length**5 \
           + cutoff_dxx(x_length)*(-18*x**2*y**2 + 3*y**4)/x_length**6 \
           + cutoff_dx(x_length)*(3*x**4 + 18*x**2*y**2)/x_length**7

def dxxx_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dxxx(x_length)*x**3./x_length**3 \
           + cutoff_dxx(x_length)*(2*y**2*x)/x_length**4 \
           - cutoff_dx(x_length)*(3*x*y**2)/x_length**5

def dyyy_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dxxx(x_length)*y**3/x_length**3 \
           + cutoff_dxx(x_length)*(3*x**2*y)/x_length**4 \
           - cutoff_dx(x_length)*(3*y*x**2)/x_length**5

def dxxy_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dxxx(x_length)*x**2*y/x_length**3 \
           + cutoff_dxx(x_length)*(-2*y*x**2 + y**3)/x_length**4 \
           + cutoff_dx(x_length)*(2*y*x**2 - y**3)/x_length**5

def dxyy_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dxxx(x_length)*x*y**2./x_length**3 \
           + cutoff_dxx(x_length)*(-2*x*y**2 + x**3)/x_length**4 \
           + cutoff_dx(x_length)*(2*x*y**2 - x**3)/x_length**5
    #return cutoff_dxxx(x_length)*x*y**2./x_length**3 \
    #       - cutoff_dxx(x_length)*(3*x*y**2)/x_length**4 \
    #       + cutoff_dx(x_length)*(3*x*y)/x_length**5

def dxx_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dxx(x_length)*x**2/x_length**2 \
           + cutoff_dx(x_length)*y**2/x_length**3

def dyy_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dxx(x_length)*y**2/x_length**2 \
           + cutoff_dx(x_length)*x**2/x_length**3

def dxy_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dxx(x_length)*x*y/x_length**2 \
           - cutoff_dx(x_length)*x*y/x_length**3

def dx_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dx(x_length)*x/x_length

def dy_phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff_dx(x_length)*y/x_length

def phi(x, y):
    x_length = np.sqrt(x**2 + y**2)
    if isinstance(x_length, (float, int)):
        x_length = max(0.5, min(abs(x_length), 7./8))
    else:
        n = x_length.shape
        x_length = np.maximum(0.5*np.ones(n), np.minimum(np.absolute(x_length), 7./8*np.ones(n)))
    return cutoff(x_length)




def dyyyy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta2**4)*N

def dxxxx_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta1**4)*N

def dxxyy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta1**2*eta2**2)*N

def dxxxy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta1**3*eta2)*N

def dxyyy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta1*eta2**3)*N

def dtxxy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3*eta1**2*eta2)*N

def dtxyy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3*eta1*eta2**2)*N

def dtyyy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3*eta2**3)*N

def dtxxx_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3*eta1**3)*N

def dttxy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3**2*eta1*eta2)*N

def dttxx_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3**2*eta1**2)*N

def dttyy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3**2*eta2**2)*N

def dtxx_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3*eta1**2)

def dtyy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3*eta2**2)

def dtxy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3*eta1*eta2)

def dxxx_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta1**3)

def dyyy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta2**3)

def dxxy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta1**2*eta2)

def dxyy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta1*eta2**2)

def dttx_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3**2*eta1)

def dtty_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3**2*eta2)

def dxx_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta1**2)/N

def dyy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta2**2)/N

def dxy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta1*eta2)/N

def dtx_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3*eta1)/N

def dty_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3*eta2)/N

def dx_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta1)/N**2

def dy_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta2)/N**2

def fpsi_b(N, eta1, eta2, eta3 ):
    return 1./N**3

def dyyt_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta2**2*eta3)

def dxyt_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta1*eta2*eta3)

def dxxt_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta1**2*eta3)

def dxt_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta1*eta3)/N

def dyt_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta2*eta3)/N

def dt_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3)/N**2

def dtt_psi_b(N, eta1, eta2, eta3 ):
    return np.abs(eta3**2)/N

def dxxyy_phi_b():
    return 31480

def dxxxy_phi_b():
    return 31480

def dxyyy_phi_b():
    return 31480

def dyyyy_phi_b():
    return 31480

def dxxxx_phi_b():
    return 31480

def dxxx_phi_b():
    return 2004

def dyyy_phi_b():
    return 2004

def dxxy_phi_b():
    return 2004

def dxyy_phi_b():
    return 2004

def dxx_phi_b():
    return 84

def dyy_phi_b():
    return 84

def dxy_phi_b():
    return 84

def dx_phi_b():
    return 7

def dy_phi_b():
    return 7

def phi_b():
    return 1

def cutoff_b():
    return 1

def cutoff_dx_b():
    return 7

def cutoff_dxx_b():
    return 84

def cutoff_dxxx_b():
    return 2004

def cutoff_dxxxx_b():
    return 31480
