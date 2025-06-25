function [symbols, shift_comps, res] = shift_observer2_decider(sequence, SF, B, OSR, LDRO, initial_shift, initial_rate)
    % initial shift in Hz, initial rate in Hz/second

    if LDRO==true
        kp = 0.4;
        ki = 0.04 ;
        ZP = 2;
    else
        kp = 0.1;
        ki = 0.01;
        ZP = 8;
    end

    symbol_len = (2^SF*OSR);
    n_symbols = length(sequence)/symbol_len;
    symbols = zeros(1, n_symbols);
    dc = downchirp(SF, B, 1);
    
    if LDRO==false
        shift_comp = initial_shift*2^SF/B;
    else
        shift_comp = 0.25*initial_shift*2^SF/B;
    end
    shift_comps = zeros(1, n_symbols);
    res = zeros(1, n_symbols);
    rates = zeros(1, n_symbols);
        
    y_int1_1 = 0;

    y_int2_1 = shift_comp;

    for i=1:n_symbols
        symbol_ds = downsample(sequence((i-1)*symbol_len+1:i*symbol_len), OSR);
        dechirped = symbol_ds.*dc;
        dechirped_zp = [dechirped zeros(1, (ZP - 1)*length(dechirped))];

        if LDRO==false
            fftres = circshift(abs(fft(dechirped_zp)), myround(-shift_comp*ZP));
        else
            fftres = circshift(abs(fft(dechirped_zp)), myround(-shift_comp*4*ZP));
        end

        [~, maxind] = max(downsample(fftres, ZP));
        s = maxind - 1;
        
        [~, maxind_zp] = max(fftres);
        s_zp = (maxind_zp - 1)/ZP;

        if LDRO == true
            symbols(i) = mod(myround(s/4), 2^(SF-2));
            x_int1 = (s_zp/4 - myround(s_zp/4));
        else
            symbols(i) = mod(s, 2^SF);
            x_int1 = (s_zp - myround(s_zp));
        end

        y_int1 = ki*x_int1 + y_int1_1;
        y_int1_1 = y_int1;
        
        x_int2 = y_int1 + kp*x_int1;
        y_int2 = x_int2 + y_int2_1;
        y_int2_1 = y_int2;
    
        rates(i) = y_int1;
        res(i) = x_int1;
        shift_comps(i) = shift_comp;
        shift_comp = y_int2;
    end
end