function [errors, total_symbols] = sat_test_param(decider, ELEV, SNR, SF, LDRO, CENTER_FREQ)
    B = 125000;
    GB = 128000;
    OSR = 2;
    N_SYMBOLS = 50;
    Ts = 2^SF/B;
    T = N_SYMBOLS*Ts;
    ALTITUDE = 500e3;
    N_SAMPLES = N_SYMBOLS*2^SF*OSR;
    N_RUNS = 2000;
    N_PAR = 100;
    N_ERR = 50;
    CONV = 0.01;
    MIN_SER = 1e-3;

    Hd_filter = filter_design(B, OSR, GB);
    Hd = Hd_filter.Numerator;
    filter_delay = floor(length(Hd)/2);

    basic_diffs = zeros(N_PAR, N_SYMBOLS);

    fprintf("running elevation %f deg; snr %f\n", ELEV*180/pi, SNR);
    shifts = satellite_shifts(ALTITUDE, CENTER_FREQ, ELEV, B*OSR, N_SAMPLES);
    initial_shift = shifts(1);
    initial_rate = (shifts(11) - shifts(1))*(10*B*OSR);
    
    errors = 0;
    total_symbols = 0;
    for j=1:N_RUNS
        fprintf("run %d/%d\n", j, N_RUNS)
        parfor k=1:N_PAR
            [sequence, symbols] = generate_symbol_sequence(N_SYMBOLS, SF, B, OSR, LDRO);
            noise = (sqrt(OSR)/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20) + (sqrt(OSR)*1j/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20);
            sequence_shifted = shiftfun(sequence, shifts, B*OSR) + noise;
            flt_out = filter(Hd, 1, [sequence_shifted zeros(1, filter_delay)]);
            sequence_flt = flt_out(filter_delay+1:end);
            out = decider(sequence_flt, SF, B, OSR, LDRO, initial_shift, initial_rate);
            basic_diffs(k, :) = abs(symbols - out);
        end
        errors = errors + sum(sign(results(basic_diffs, SF, LDRO)));
        total_symbols = total_symbols + N_PAR*N_SYMBOLS;
    
        [ser, ci] = binofit((errors), total_symbols);
        ci_rw_log = abs((log(ci(2)) - log(ci(1)))/log(MIN_SER));
        fprintf("errors = %d; ser = %f; ci w = %f; ci_rw_log = %f\n", (errors), ser, ci(2) - ci(1), ci_rw_log);
        if (errors) > N_ERR && ci_rw_log < CONV
            break;
        end
        if total_symbols > 10/MIN_SER && (ci(2) < MIN_SER || (errors) == 0 || (errors) > 200000)
            fprintf("SER CI upper bound lower than %f, skipping\n", MIN_SER);
            break;
        end
    end
end