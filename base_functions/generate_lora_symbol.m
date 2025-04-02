function samples = generate_lora_symbol(s, SF, B, OSR, shift)
    M = 2^SF;
    Ts = M/B;
    Th = Ts - s*Ts/M;
    t = linspace(0, Ts - Ts/(M*OSR), M*OSR);

    wrap = OSR*(M-s+1);
    phases = (t.*(-B/2 + s*B/M) + (t.^2).*B/(2*Ts)) + shift*t;
    phases(wrap:end) = phases(wrap:end) - (t(wrap:end) - Th).*(B);
    samples = exp(2.0j*pi*phases);
end