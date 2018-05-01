function [y] = dx_phi(x, y)
    x_length = sqrt(x.^2 + y.^2);
    x_length = max(0.5,min(abs(x_length),7/8));
    y = cutoff_dx(x_length).*x./x_length;