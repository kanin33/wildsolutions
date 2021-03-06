function [y] = dxxx_phi(x, y)
    x_length = sqrt(x.^2 + y.^2);
    x_length = max(0.5,min(abs(x_length),7/8));
    y = cutoff_dxxx(x_length).*x.^3./x_length.^3 ...
           + cutoff_dxx(x_length).*(2.*x.^2.*y)./x_length.^4 ...
           - cutoff_dx(x_length).*(3.*x.*y.^2)./x_length.^5;