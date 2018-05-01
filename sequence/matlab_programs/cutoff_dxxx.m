function [y] = cutoff_dxxx(x)
    x = max(0.5,min(abs(x),7/8));
    y = -524*(16*x - 11).*(96.*x - 66).*(16*x.^2 - 22*x + 7).^2 ...
        - 8384*(16*x.^2 - 22*x + 7).^3;