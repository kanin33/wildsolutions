function [y] = cutoff_dx(x)
    x = max(0.5,min(abs(x),7/8));
    y = -(143360/2187)*(16*x.^2 - 22*x + 7).^4;