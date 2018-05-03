import warnings
warnings.filterwarnings("ignore")
import numpy as np

"""
This module computes the constant subsolution (v_tilde, u_tilde)
(When the standard solution consists of a 1-shock and a 3-shock)
"""

def find_mu(rho_m, rho_p, rho_1, v_m, v_p, p):
    """
    Return mu
    """

    mu = (rho_m*v_m[1] - rho_p*v_p[1])/float(rho_m - rho_p)
    tmp_sqrt = rho_p*rho_m*(v_m[1] - v_p[1])**2 \
               - (rho_m - rho_p)*(p(rho_m) - p(rho_p))
    mu_0 = mu - 1/float(rho_m - rho_p) \
                *np.sqrt(tmp_sqrt*(rho_1 - rho_p)/float(rho_1 - rho_m))
    mu_1 = mu - 1/float(rho_m - rho_p) \
                *np.sqrt(tmp_sqrt*(rho_1 - rho_m)/float(rho_1 - rho_p))

    return mu_0, mu_1

def find_delta_1(rho_m, rho_p, rho_1, v_m, v_p, p):
    """
    Return delta_1
    """

    tmp_sqrt = rho_p*rho_m*(v_m[1] - v_p[1])**2 \
               - (rho_m - rho_p)*(p(rho_m) - p(rho_p))
    tmp_sqrt *= (rho_1 - rho_p)/float(rho_1 - rho_m)
    tmp_sqrt = np.sqrt(tmp_sqrt)
    delta_1 = -(p(rho_1) - p(rho_m))/float(rho_1) \
              + rho_m*float(rho_1 - rho_m)/(rho_1**2*(rho_m - rho_p)**2)\
              *(rho_p*(v_p[1] - v_m[1]) + tmp_sqrt)**2
    return delta_1

def find_v12(rho_m, rho_p, rho_1, v_m, v_p, p):
    """
    Return v_12
    """

    tmp_sqrt = rho_p*rho_m*(v_m[1] - v_p[1])**2 \
               - (rho_m - rho_p)*(p(rho_m) - p(rho_p))
    tmp_sqrt *= (rho_1 - rho_m)*(rho_1 - rho_p)
    tmp_sqrt = np.sqrt(tmp_sqrt)

    v_12 = 1./(rho_1*(rho_m - rho_p))*(rho_m*v_m[1]*(rho_1 - rho_p)
            - rho_p*v_p[1]*(rho_1 - rho_m) - tmp_sqrt)
    return v_12

def find_v_u_tilde_C(v_m, v_12, delta):
    """
    Return v_tilde, u_tilde and C_1
    """
    delta_1 = delta[0]; delta_2 = delta[1]
    v_11 = v_m[0]
    u_112 = v_11*v_12
    C_1 = delta_1 + delta_2 + v_11**2 + v_12**2
    u_111 = (delta_2 - delta_1 + v_11**2 - v_12**2)/2.

    v_tilde = np.array([v_11, v_12])
    u_tilde = np.array([[u_111, u_112],
               [-u_112, -u_111]])

    return v_tilde, u_tilde, C_1

