function [y] = dyyy_phi(x, y)
    x_length = sqrt(x.^2 + y.^2);
    x_length = max(0.5,min(abs(x_length),7/8));
    y = cutoff_dxxx(x_length).*y.^3./x_length.^3 ...
           + cutoff_dxx(x_length).*(2.*y.^2.*x)./x_length.^4 ...
           - cutoff_dx(x_length).*(3.*y.*x.^2)./x_length.^5;