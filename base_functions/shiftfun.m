function [y, ps] = shiftfun(x, fqs, fs)
    lx = length(x);
    T = lx/fs;
    t = linspace(0, T - T/lx, lx);
    ps = 2.0*pi*cumtrapz(fqs)/fs;
    y = x.*exp(1j*ps);
end