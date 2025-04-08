function sat_test(SF, LDRO, CENTER_FREQ)
    clf;
    % SF = 12;
    % LDRO = true;
    B = 125000;
    GB = 128000;
    OSR = 2;
    SNR = -9;
    N_SYMBOLS = 50;
    Ts = 2^SF/B;
    T = N_SYMBOLS*Ts;
    ALTITUDE = 500e3;
    % CENTER_FREQ = 915e6;
    N_SAMPLES = N_SYMBOLS*2^SF*OSR;
    ELEVATION = -5*pi/180;
    N_RUNS = 400;
    N_PAR = 100;
    MIN_ERR = 2000;
    
    elevs = [0, 45*pi/180, 89.9*pi/180];
    snrs = linspace(-12, -9, 5) - 3*(SF - 7);
    % elevs = [89.8*pi/180,];
    % snrs = [-12, -10.5, -9] - 3*(SF - 7);
    configs = combinations(elevs, snrs);
    basic_bers = zeros(1, height(configs));
    basic_sers = zeros(1, height(configs));
    observer_bers = zeros(1, height(configs));
    observer_sers = zeros(1, height(configs));
    
    for j=1:height(configs)
        fprintf("config %d/%d\n", j, height(configs));
        elev = configs(j, :).elevs;
        snr = configs(j, :).snrs;
        shifts = satellite_shifts(ALTITUDE, CENTER_FREQ, elev, B*OSR, N_SAMPLES);
        initial_shift = shifts(1);
        initial_rate = (shifts(11) - shifts(1))*(10*B*OSR);
        
        Hd_filter = filter_design(B, OSR, GB);
        Hd = Hd_filter.Numerator;
        filter_delay = floor(length(Hd)/2);
        
        basic_diffs = zeros(N_PAR, N_RUNS, N_SYMBOLS);
        observer_diffs = zeros(N_PAR, N_RUNS, N_SYMBOLS);
        
        [sequence, symbols] = generate_symbol_sequence(N_SYMBOLS, SF, B, OSR, LDRO);
        total_runs = 0;
        for k = 1:N_RUNS
            parfor i = 1:N_PAR
                noise = (1/sqrt(2))*randn(1, N_SAMPLES)*10^(-snr/20) + (1j/sqrt(2))*randn(1, N_SAMPLES)*10^(-snr/20);
                sequence_shifted = shift(sequence, shifts, B*OSR) + noise;
                flt_out = filter(Hd, 1, [sequence_shifted zeros(1, filter_delay)]);
                sequence_shifted_flt = flt_out(filter_delay+1:end);
                % fprintf("%d\n", symbols(1));
                basic_out = basic_decider(sequence_shifted_flt, SF, B, OSR, LDRO, initial_shift, initial_rate);
                observer_out = shift_observer2_decider(sequence_shifted_flt, SF, B, OSR, LDRO, initial_shift, initial_rate);
                basic_diffs(i, k, :) = abs(symbols - basic_out);
                observer_diffs(i, k, :) = abs(symbols - observer_out);
            end
            fprintf("run %d/%d\n", k, N_RUNS);
    
            basic_res = results(basic_diffs, SF, LDRO);
            observer_res = results(observer_diffs, SF, LDRO);
            total_runs = total_runs + length(symbols)*N_PAR;
            
            if sum(sign(basic_res)) > MIN_ERR && sum(sign(observer_res)) > MIN_ERR
                fprintf("skipping next runs as error count is > %d for both\n", MIN_ERR);
                break;
            end
        end
        
        res = results(basic_diffs, SF, LDRO);
        ber = sum(res)/(total_runs*SF);
        ser = sum(sign(res))/(total_runs);
        fprintf("basic ber = %f%%; ser = %f%%\n", ber*100, ser*100);
        basic_sers(j) = ser;
        basic_bers(j) = ber;
    
        res = results(observer_diffs, SF, LDRO);
        ber = sum(res)/(total_runs*SF);
        ser = sum(sign(res))/(total_runs);
        fprintf("observer ber = %f%%; ser = %f%%\n", ber*100, ser*100);
        observer_sers(j) = ser;
        observer_bers(j) = ber;
    end
    
    
    clf;
    hold on;
    ylim([1e-4, 1]);
    
    for i=1:length(elevs)
        basic_sers_r = zeros(1, length(snrs));
        snrs_r = zeros(1, length(snrs));
        k = 1;
        for j=1:height(configs)
            if configs(j, :).elevs == elevs(i)
                snrs_r(k) = configs(j, :).snrs;
                basic_sers_r(k) = basic_sers(j);
                k = k + 1;
            end
        end
        t = sprintf("Basic - Initial Elev. = %.1f °", elevs(i)*180/pi);
        semilogy(snrs, basic_sers_r, "DisplayName", t);
    end
    
    for i=1:length(elevs)
        observer_sers_r = zeros(1, length(snrs));
        snrs_r = zeros(1, length(snrs));
        k = 1;
        for j=1:height(configs)
            if configs(j, :).elevs == elevs(i)
                snrs_r(k) = configs(j, :).snrs;
                observer_sers_r(k) = observer_sers(j);
                k = k + 1;
            end
        end
        t = sprintf("Observer - Initial Elev. = %.1f °", elevs(i)*180/pi);
        semilogy(snrs, observer_sers_r, "DisplayName", t);
    end
    
    set(gca, 'YScale', 'log');
    legend(gca, "Location", "best");
    grid();
    ylabel("SER");
    xlabel("SNR");
    t = sprintf("Satellite Altitude = %.0f km; Center Frequency = %.0f MHz; \nSF=%d; Symbols/Packet=%d; BW = %.0f kHz", ALTITUDE/1e3, CENTER_FREQ/1e6, SF, N_SYMBOLS, B/1e3);
    title(t);
    if LDRO == true
        ldro_str = "ldro";
    else
        ldro_str = "";
    end
    t = sprintf("sat_test_%d_%d%s.png", round(CENTER_FREQ/1e6), SF, ldro_str);
    saveas(gcf(), t, "png");
    hold off;
    % save;
end