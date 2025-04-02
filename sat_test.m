SF = 10;
B = 125000;
GB = 128000;
OSR = 2;
SNR = -9;
N_SYMBOLS = 50;
Ts = 2^SF/B;
T = N_SYMBOLS*Ts;
ALTITUDE = 500e3;
CENTER_FREQ = 915e6;
N_SAMPLES = N_SYMBOLS*2^SF*OSR;
ELEVATION = -5*pi/180;
N_RUNS = 10;
N_PAR = 1000;

n_snrs = 5;n_elevs = 3;
elevs = [ -20*pi/180, -20*pi/180, -20*pi/180, -20*pi/180, -20*pi/180, ...
           -5*pi/180,  -5*pi/180,  -5*pi/180,  -5*pi/180,  -5*pi/180, ...
                   0,          0,          0,          0,          0];

snrs = [-15, -13.5, -12, -10.5, -9, -15, -13.5, -12, -10.5, -9, -15, -13.5, -12, -10.5, -9] - 3*(SF - 7);
basic_bers = zeros(1, length(elevs));
basic_sers = zeros(1, length(elevs));
observer_bers = zeros(1, length(elevs));
observer_sers = zeros(1, length(elevs));

for j=1:length(elevs)
    fprintf("config %d/%d\n", j, length(elevs));
    elev = elevs(j);
    snr = snrs(j);
    shifts = satellite_shifts(ALTITUDE, CENTER_FREQ, elev, B*OSR, N_SAMPLES);
    initial_shift = shifts(1);
    initial_rate = (shifts(11) - shifts(1))*(10*B*OSR);
    
    Hd_filter = filter_design(B, OSR, GB);
    Hd = Hd_filter.Numerator;
    filter_delay = floor(length(Hd)/2);
    
    basic_diffs = zeros(N_PAR, N_RUNS, N_SYMBOLS);
    observer_diffs = zeros(N_PAR, N_RUNS, N_SYMBOLS);
    
    [sequence, symbols] = generate_symbol_sequence(N_SYMBOLS, SF, B, OSR);
    
    for k = 1:N_RUNS
        parfor i = 1:N_PAR
            noise = (1/sqrt(2))*randn(1, N_SAMPLES)*10^(-snr/20) + (1j/sqrt(2))*randn(1, N_SAMPLES)*10^(-snr/20);
            sequence_shifted = shift(sequence, shifts, B*OSR) + noise;
            flt_out = filter(Hd, 1, [sequence_shifted zeros(1, filter_delay)]);
            sequence_shifted_flt = flt_out(filter_delay+1:end);
            
            basic_out = basic_decider(sequence_shifted_flt, SF, B, OSR, initial_shift, initial_rate);
            observer_out = shift_observer_decider(sequence_shifted_flt, SF, B, OSR, initial_shift, initial_rate);
            basic_diffs(i, k, :) = abs(symbols - basic_out);
            observer_diffs(i, k, :) = abs(symbols - observer_out);
        end
        fprintf("run %d/%d\n", k, N_RUNS);
    end
    
    res = results(basic_diffs, SF);
    ber = sum(res)/(length(res)*SF);
    ser = sum(sign(res))/(length(res));
    fprintf("basic ber = %f%%; ser = %f%%\n", ber*100, ser*100);
    basic_sers(j) = ser;
    basic_bers(j) = ber;

    res = results(observer_diffs, SF);
    ber = sum(res)/(length(res)*SF);
    ser = sum(sign(res))/(length(res));
    fprintf("observer ber = %f%%; ser = %f%%\n", ber*100, ser*100);
    observer_sers(j) = ser;
    observer_bers(j) = ber;
end


clf;
hold on;
ylim([1e-3, 1]);
basic_sers = reshape(basic_sers, n_snrs, n_elevs)';
basic_bers = reshape(basic_bers, n_snrs, n_elevs)';
snrs = reshape(snrs, n_snrs, n_elevs)';
elevs = reshape(elevs, n_snrs, n_elevs)';
for i=1:n_elevs
    t = sprintf("Basic - Initial Elev. = %.1f °", elevs(i, 1)*180/pi);
    semilogy(snrs(i,:), basic_sers(i, :), "DisplayName", t);
end

observer_sers = reshape(observer_sers, n_snrs, n_elevs)';
observer_bers = reshape(observer_bers, n_snrs, n_elevs)';
for i=1:n_elevs
    t = sprintf("Observer - Initial Elev. = %.1f °", elevs(i, 1)*180/pi);
    semilogy(snrs(i,:), observer_sers(i, :), "DisplayName", t);
end
set(gca, 'YScale', 'log');
legend();
grid();
ylabel("SER");
xlabel("SNR");
t = sprintf("Satellite Altitude = %.0f km; Center Frequency = %.0f MHz; \nSF=%d; Symbols/Packet=%d; BW = %.0f kHz", ALTITUDE/1e3, CENTER_FREQ/1e6, SF, N_SYMBOLS, B/1e3);
title(t);
saveas(gcf(), "sat_test.png", "png");
hold off;