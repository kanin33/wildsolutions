function [z] = dy_phi(x, y)
    x_length = sqrt(x.^2 + y.^2);
    x_length = max(0.5,min(abs(x_length),7/8));
    z = cutoff_dx(x_length).*y./x_length;