function symbols = shift_observer2_decider(sequence, SF, B, OSR, LDRO, initial_shift, initial_rate)
    % initial shift in Hz, initial rate in Hz/second
    kp = 0.3; ki = 0.1; c = 1;

    ZP = 32;
    symbol_len = (2^SF*OSR);
    n_symbols = length(sequence)/symbol_len;
    symbols = zeros(1, n_symbols);
    dc = downchirp(SF, B, 1);
    
    shift_comp = -initial_shift*2^SF/B;
        
    y_int1_1 = 0;

    y_int2_1 = shift_comp;

    for i=1:n_symbols
        symbol_ds = downsample(sequence((i-1)*symbol_len+1:i*symbol_len), OSR);
        dechirped = symbol_ds.*dc;
        dechirped_zp = [dechirped zeros(1, (ZP - 1)*length(dechirped))];
        fftres = abs(fft(dechirped_zp));

        [~, maxind] = max(downsample(circshift(fftres, round(shift_comp*ZP)), ZP));
        s = maxind - 1;

        if LDRO == true
            symbols(i) = mod(round(s/4), 2^(SF-2));
        else
            symbols(i) = mod(round(s), 2^SF);
        end

        [~, maxind_zp] = max(fftres);
        s_zp = (maxind_zp - 1)/ZP + shift_comp;
        
        x_int1 = round(s_zp) - s_zp;
        y_int1 = x_int1 + y_int1_1;
        y_int1_1 = y_int1;
        
        x_int2 = y_int1*ki + kp*x_int1;
        y_int2 = c*x_int2 + y_int2_1;
        y_int2_1 = y_int2;

        shift_comp = y_int2;
    end
end