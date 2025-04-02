function y = shift(x, freqs, fs)
    lx = length(x);
    T = lx/fs;
    t = linspace(0, T - T/lx, lx);
    y = x.*exp(2.0j*pi.*freqs.*t);
end