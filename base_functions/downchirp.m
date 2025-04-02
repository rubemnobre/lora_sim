function samples = downchirp(SF, B, OSR)
    M = 2^SF;
    Ts = M/B;
    t = linspace(0, Ts - Ts/(M*OSR), M*OSR);
    phases = (t.*(B/2) - (t.^2).*B/(2*Ts));
    samples = exp(2.0j*pi*phases);
end