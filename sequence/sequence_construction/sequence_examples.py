from gen_sequence import *
from subsolution import *
from plot_functions import *

def example_5_6():
    # Example 5.6
    rho_m   = 1; rho_p = 4; rho_1 = 8
    v_m     = [0, 8./3.*np.sqrt(10)]; v_p = [0, -5./6*np.sqrt(13)]
    p       = lambda rho: rho**2
    delta_2 = 1
    mu      = find_mu(rho_m, rho_p, rho_1, v_m, v_p, p)
    delta   = [find_delta_1(rho_m, rho_p, rho_1, v_m, v_p, p), delta_2]
    v_12    = find_v12(rho_m, rho_p, rho_1, v_m, v_p, p)
    v_tilde, u_tilde, C_1 = find_v_u_tilde_C(v_m, v_12, delta)

    return v_m, v_p, mu, delta, v_12, v_tilde, u_tilde, C_1

def print_ex_5_6_sub(plot_sub=True):
    """
    Plot subsolution
    """

    v_m, v_p, mu, delta, v_12, v_tilde, u_tilde, C_1 = example_5_6()
    print "mu_0: %.2f\nmu_1: %.2f" % (mu[0], mu[1])
    print "delta_1: %.2f\ndelta_2: %.2f" % (delta[0], delta[1])
    print "v_12: %.2f" % v_12
    print "v_tilde: [%.2f, %.2f]" % (v_tilde[0], v_tilde[1])
    print "u_tilde:"
    print "[[%.2f, %.2f]]" % (u_tilde[0][0], u_tilde[0][1])
    print "[[%.2f, %.2f]]" % (u_tilde[1][0], u_tilde[1][1])
    print "C_1: %.2f" % C_1

def compute_gen_ex_5_6(max_iterations, t_min, t_max,
                                       x_1_min, x_1_max,
                                       x_2_min=-1, x_2_max=1,
                                       plot_v=False, t=1.5,
                                       dist_filename="dist.txt",
                                       r_values=None,
                                       dirname="test_dir",
                                       figname="test.png",
                                       quiver=False,
                                       stream=False,
                                       num_levels=30,
                                       contour=False,
                                       curl=False,
                                       stream_plt=False,
                                       show=True,
                                       dpi=None,
                                       check_steps=True,
                                       z_k=None,
                                       statistics_list=None
                       ):

    """
    Compute a sequence corresponding to values in example 5.6
    """
    t = np.float64(t)
    seq, z_k = generate_sequence(example_5_6, t_min=t_min, t_max=t_max,
                                              x_1_min=x_1_min, x_1_max=x_1_max,
                                              max_iterations=max_iterations,
                                              dist_filename=dist_filename,
                                              r_values=r_values,
                                              dirname=dirname,
                                              check_steps=check_steps,
                                              z_k=z_k,
                                              statistics_list=statistics_list)
    v = lambda t, x1, x2: z_k(t, x1, x2)[0]

    if plot_v:
        plot_v_x1_x2(v, x_2_min, x_2_max,
                        x_1_min, x_1_max,
                        t, n=1000,
                        quiver=quiver,
                        stream=stream,
                        num_levels=num_levels,
                        contour=contour,
                        curl=curl,
                        stream_plt=stream_plt,
                        show=show,
                        filename=figname,
                        dpi=dpi
                     )
    return seq, z_k
