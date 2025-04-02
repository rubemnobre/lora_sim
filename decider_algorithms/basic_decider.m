function symbols = basic_decider(sequence, SF, B, OSR, initial_shift, initial_rate)
    % initial shift in Hz, initial rate in Hz/second
    symbol_len = (2^SF*OSR);
    n_symbols = length(sequence)/symbol_len;
    symbols = zeros(1, n_symbols);
    dc = downchirp(SF, B, 1);
    
    Ts = 2^SF/B;
    t = linspace(0, Ts - Ts/2^SF, 2^SF);
    int_shift = round(initial_shift/B);
    frac_shift = initial_shift/B - int_shift;
    
    frac_scomp = exp(-2j*pi*B*frac_shift.*t);
    
    for i=1:n_symbols
        symbol_ds = downsample(sequence((i-1)*symbol_len+1:i*symbol_len), OSR);
        dechirped = symbol_ds.*frac_scomp.*dc;
        fftres = abs(fft(dechirped));
        [m, maxind] = max(fftres);
        symbols(i) = mod(maxind - 1 - int_shift, 2^SF);
    end
end