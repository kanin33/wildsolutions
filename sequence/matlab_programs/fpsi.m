function [y] = fpsi(N, eta1, eta2, eta3, x1_new, x2_new, t_new)
    y = -sin(N*(eta1*x1_new + eta2*x2_new + eta3*t_new))/N^3;