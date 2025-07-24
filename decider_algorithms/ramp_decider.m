function [symbols, shift_comps] = ramp_decider(sequence, SF, B, OSR, LDRO, initial_shift, preamble_present)
    % initial shift in Hz, initial rate in Hz/second
    SW = 52;
    N = 10;
    if LDRO==true
        ZP = 2;
    else
        ZP = 8;
    end

    symbol_len = (2^SF*OSR);
    if preamble_present
        n_symbols = round(length(sequence)/symbol_len - 12.25);
    else
        n_symbols = length(sequence)/symbol_len;
    end
    symbols = zeros(1, n_symbols);
    dc = downchirp(SF, B, 1);
    
    if LDRO==false
        shift_comp = initial_shift*2^SF/B;
    else
        shift_comp = 0.25*initial_shift*2^SF/B;
    end
    shift_comps = zeros(1, n_symbols);

    offset = 0;

    % Run loop on the 8 upchirps of the preamble
    if preamble_present
        offset = round(symbol_len*12.25);
        expected = zeros(1, 10);
        s1 = bitshift(SW, -4)*8;
        s2 = bitand(SW, 0xF)*8;
        expected(9) = s1;
        expected(10) = s2;
        vals = zeros(1, 10);
        for i=1:10
            symbol_ds = downsample(sequence((i-1)*symbol_len+1:i*symbol_len), OSR);
            dechirped = symbol_ds.*dc;
            dechirped_zp = [dechirped zeros(1, (ZP - 1)*length(dechirped))];
    
            if LDRO==false
                sc = myround(-shift_comp*ZP);
            else
                sc = myround(-shift_comp*4*ZP);
            end
            fftres = abs(fft(dechirped_zp));
    
            [~, maxind_zp] = max(fftres);
            s_zp = (maxind_zp + sc - 1)/ZP;
            
            if LDRO == true
                s = (s_zp/4 - expected(i)/4);
                if s > 2^SF/8
                    s = s - 2^SF/4;
                end
            else
                s = (s_zp - expected(i));
                if s > 2^SF/2
                    s = s - 2^SF;
                end
            end
            vals(i) = s;
        end
        last_val = vals(1);
        last_idx = 1;
        idx = 2;
        rates = zeros(1,0);
        while idx < 10
            if vals(idx) <= 4 && vals(idx) >= -4
                rates = [rates (vals(idx) - last_val)/(idx - last_idx)];
                last_val = vals(idx);
                last_idx = idx;
            end
            idx = idx + 1;
        end
        if isempty(rates)
            rate_est = 0;
        else
            rate_est = mean(rates);
        end
        shift_comp = shift_comp + rate_est*12.25;
    end
    
    vals = zeros(1, N);
    vi = 1;
    for i=1:n_symbols
        symbol_ds = downsample(sequence((i-1)*symbol_len+1+offset:i*symbol_len+offset), OSR);
        dechirped = symbol_ds.*dc;
        dechirped_zp = [dechirped zeros(1, (ZP - 1)*length(dechirped))];

        if LDRO==false
            sc = myround(-shift_comp*ZP);
        else
            sc = myround(-shift_comp*4*ZP);
        end

        fftres = abs(fft(dechirped_zp));

        [~, maxind_zp] = max(fftres);
        s_zp = (maxind_zp + sc - 1)/ZP;

        if LDRO == true
            symbols(i) = mod(myround(s_zp/4), 2^(SF-2));
            vals(vi) = (s_zp/4) - myround(s_zp/4);
        else
            symbols(i) = mod(myround(s_zp), 2^SF);
            vals(vi) = (s_zp) - myround(s_zp);
        end

        if vi == N
            delta_rate = mean(diff(vals));
            rate_est = rate_est + delta_rate;
            vi = 1;
        end
   
        shift_comps(i) = shift_comp;
        shift_comp = shift_comp + rate_est;
        vi = vi + 1;
    end
end