function sat_test_elev(SNR, SF, LDRO, CENTER_FREQ)
    clf;
    B = 125000;
    GB = 128000;
    OSR = 2;
    N_SYMBOLS = 50;
    Ts = 2^SF/B;
    T = N_SYMBOLS*Ts;
    ALTITUDE = 500e3;
    N_SAMPLES = N_SYMBOLS*2^SF*OSR;
    N_RUNS = 500;
    N_PAR = 100;
    N_ERR = 100;
    CONV = 0.02;
    
    elevs = linspace(0, pi/2, 31);
    basic_sers = zeros(1, length(elevs));
    observer_sers = zeros(1, length(elevs));

    Hd_filter = filter_design(B, OSR, GB);
    Hd = Hd_filter.Numerator;
    filter_delay = floor(length(Hd)/2);

    basic_diffs = zeros(N_PAR, N_SYMBOLS);
    observer_diffs = zeros(N_PAR, N_SYMBOLS);

    for i=1:length(elevs)
        fprintf("running elevation %f deg\n", elevs(i)*180/pi);
        shifts = satellite_shifts(ALTITUDE, CENTER_FREQ, elevs(i), B*OSR, N_SAMPLES);
        initial_shift = shifts(1);
        initial_rate = (shifts(11) - shifts(1))*(10*B*OSR);
        
        errors_basic = 0;
        errors_observer = 0;
        total_symbols_basic = 0;
        total_symbols_observer = 0;
        tic;
        prev_ser = 0;
        for j=1:N_RUNS
            fprintf("basic run %d/%d\n", j, N_RUNS)
            parfor k=1:N_PAR
                rng(k + j*N_PAR);
                [sequence, symbols] = generate_symbol_sequence(N_SYMBOLS, SF, B, OSR, LDRO);
                noise = (1/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20) + (1j/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20);
                sequence_shifted = shift(sequence, shifts, B*OSR) + noise;
                flt_out = filter(Hd, 1, [sequence_shifted zeros(1, filter_delay)]);
                sequence_flt = flt_out(filter_delay+1:end);
                basic_out = basic_decider(sequence_flt, SF, B, OSR, LDRO, initial_shift, initial_rate);
                basic_diffs(k, :) = abs(symbols - basic_out);
            end
            errors_basic = errors_basic + sum(sign(results(basic_diffs, SF, LDRO)));
            total_symbols_basic = total_symbols_basic + N_PAR*N_SYMBOLS;


    
            [ser, ci] = binofit(errors_basic, total_symbols_basic);
            ci_rw_log = abs((log(ci(2)) - log(ci(1)))/log(ser));
            fprintf("errors basic = %d; ser = %f; ci w = %f; ci_rw_log = %f\n", errors_basic, ser, ci(2) - ci(1), ci_rw_log);
            if errors_basic > N_ERR && ci_rw_log < CONV
                break;
            end
            prev_ser = ser;
        end

        for j=1:N_RUNS
            fprintf("observer run %d/%d\n", j, N_RUNS);
            parfor k=1:N_PAR
                rng(k + j*N_PAR);
                [sequence, symbols] = generate_symbol_sequence(N_SYMBOLS, SF, B, OSR, LDRO);
                noise = (1/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20) + (1j/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20);
                sequence_shifted = shift(sequence, shifts, B*OSR) + noise;
                flt_out = filter(Hd, 1, [sequence_shifted zeros(1, filter_delay)]);
                sequence_flt = flt_out(filter_delay+1:end);
                observer_out = shift_observer2_decider(sequence_flt, SF, B, OSR, LDRO, initial_shift, initial_rate);
                observer_diffs(k, :) = abs(symbols - observer_out);
            end
            errors_observer = errors_observer + sum(sign(results(observer_diffs, SF, LDRO)));
            total_symbols_observer = total_symbols_observer + N_PAR*N_SYMBOLS;

            [ser, ci] = binofit(errors_observer, total_symbols_observer);
            ci_rw_log = abs((log(ci(2)) - log(ci(1)))/log(ser));
            fprintf("errors observer = %d; ser = %f; ci w = %f; ci_rw_log = %f\n", errors_observer, ser, ci(2) - ci(1), ci_rw_log);
            if errors_observer > N_ERR && ci_rw_log < CONV
                break;
            end
            prev_ser = ser;
        end
        fprintf("elev = %f; basic ser = %.6f; observer ser = %.6f\n", elevs(i)*180/pi, errors_basic/total_symbols_basic, errors_observer/total_symbols_observer);
        toc;
        basic_sers(i) = errors_basic/total_symbols_basic;
        observer_sers(i) = errors_observer/total_symbols_observer;
    end

    hold on;
    plot(elevs*180/pi, basic_sers, "DisplayName", "Basic");
    plot(elevs*180/pi, observer_sers, "DisplayName", "Observer");
    hold off;

    ylim([1e-4, 1]);
    set(gca, 'YScale', 'log');
    legend(gca, "Location", "best");
    grid();
    ylabel("SER");
    xlabel("Initial Angle (°)");
    t = sprintf("Satellite Altitude = %.0f km; Center Frequency = %.0f MHz; \nSF = %d; Symbols/Packet=%d; BW = %.0f kHz; SNR = %d dB", ALTITUDE/1e3, CENTER_FREQ/1e6, SF, N_SYMBOLS, B/1e3, SNR);
    title(t);
    if LDRO == true
        ldro_str = "ldro";
    else
        ldro_str = "";
    end
    t = sprintf("sat_test_elev_%d_%d_%d_%s", round(CENTER_FREQ/1e6), SF, SNR, ldro_str);
    saveas(gcf(), t + ".png", "png");
    save(t + ".mat");
end