function [y] = cutoff_dxx(x)
    x = max(0.5,min(abs(x),7/8));
    y = -(8*143360/2187)*(16*x.^2 - 22*x + 7).^3.*(16*x - 11);