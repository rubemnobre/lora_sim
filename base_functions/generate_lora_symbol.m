function samples = generate_lora_symbol(s, SF, B, OSR, LDRO, shift)
    if LDRO == true
        sl = s*4;
    else
        sl = s;
    end

    M = 2^SF;
    Ts = M/B;
    Th = Ts - sl*Ts/M;
    t = linspace(0, Ts - Ts/(M*OSR), M*OSR);

    wrap = OSR*(M-sl+1);
    phases = (t.*(-B/2 + sl*B/M) + (t.^2).*B/(2*Ts)) + shift*t;
    phases(wrap:end) = phases(wrap:end) - (t(wrap:end) - Th).*(B);
    samples = exp(2.0j*pi*phases);
end