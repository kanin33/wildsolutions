function [y] = dxy_phi(x, y)
    x_length = sqrt(x.^2 + y.^2);
    x_length = max(0.5,min(abs(x_length),7/8));
    y = cutoff_dxx(x_length).*x.*y./x_length ...
           - cutoff_dx(x_length).*x.*y./x_length.^3;