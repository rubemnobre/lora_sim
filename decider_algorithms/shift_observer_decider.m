function symbols = shift_observer_decider(sequence, SF, B, OSR, LDRO, initial_shift, initial_rate)
    % initial shift in Hz, initial rate in Hz/second
    kp = 0.1;

    ZP = 4;
    symbol_len = (2^SF*OSR);
    n_symbols = length(sequence)/symbol_len;
    symbols = zeros(1, n_symbols);
    dc = downchirp(SF, B, 1);
    
    shift_comp = -initial_shift*2^SF/B;
        
    x_1 = 0;
    y_1 = shift_comp;

    for i=1:n_symbols
        symbol_ds = downsample(sequence((i-1)*symbol_len+1:i*symbol_len), OSR);
        dechirped = symbol_ds.*dc;
        dechirped_zp = [dechirped zeros(1, (ZP - 1)*length(dechirped))];
        fftres = abs(fft(dechirped_zp));
        [m, maxind] = max(fftres);
        s = (maxind - 1)/ZP + shift_comp;
        if LDRO == true
            symbols(i) = mod(round(s/4), 2^(SF-2));
        else
            symbols(i) = mod(round(s), 2^SF);
        end
        x = round(s) - s;
        y = kp*0.5*(x + x_1) + y_1;
        shift_comp = y;
        x_1 = x;
        y_1 = y;
    end
end