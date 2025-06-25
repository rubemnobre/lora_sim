
ELEV = pi/2;
SNR = -16;
SF = 10;
LDRO = false;
CENTER_FREQ = 915e6;
decider = @basic_decider;
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
CONV = 0.05;
MIN_SER = 0.1e-3;

Hd_filter = filter_design(B, OSR, GB);
Hd = Hd_filter.Numerator;
filter_delay = floor(length(Hd)/2);

basic_diffs = zeros(N_PAR, N_SYMBOLS);

fprintf("running elevation %f deg; snr %f\n", ELEV*180/pi, SNR);
shifts = satellite_shifts(ALTITUDE, CENTER_FREQ, ELEV, B*OSR, N_SAMPLES);

initial_shift = shifts(1);
initial_rate = (shifts(11) - shifts(1))*(10*B*OSR);

errors = zeros(1, N_SYMBOLS);
total_symbols = 0;
for j=1:N_RUNS
    fprintf("run %d/%d\n", j, N_RUNS)
    for k=1:N_PAR
        rng(k + j*N_PAR);
        [sequence, symbols] = generate_symbol_sequence(N_SYMBOLS, SF, B, OSR, LDRO);
        noise = (sqrt(OSR)/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20) + (sqrt(OSR)*1j/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20);
        sequence_shifted = shiftfun(sequence, shifts, B*OSR) + noise;
        flt_out = filter(Hd, 1, [sequence_shifted zeros(1, filter_delay)]);
        sequence_flt = flt_out(filter_delay+1:end);
        out = decider(sequence_flt, SF, B, OSR, LDRO, initial_shift, initial_rate);
        basic_diffs(k, :) = abs(symbols - out);
    end
    errors = errors + sum(reshape(sign(results(basic_diffs, SF, LDRO)), N_PAR, N_SYMBOLS));
    total_symbols = total_symbols + N_PAR*N_SYMBOLS;
    ci_ws = zeros(1, N_SYMBOLS);
    not_good = false;
    for k=1:N_SYMBOLS
        [ser, ci] = binofit(errors(k), total_symbols/N_SYMBOLS);
        ci_rw_log = abs(ci(2) - ci(1));
        % fprintf("errors = %d; ser = %f; ci w = %f; ci_rw_log = %f\n", errors(end), ser, ci(2) - ci(1), ci_rw_log);
        if not(ci_rw_log < CONV)
            not_good = true;
        end
        ci_ws(k) = ci_rw_log;
    end
    fprintf("max ci_ws = %f\n", max(ci_ws));
    if not_good == false
        break;
    end
end

clf();
out = errors./(total_symbols/N_SYMBOLS);
yyaxis left;
plot(downsample(shifts - shifts(1), (2^SF*OSR))./(B/(2^SF)), "Marker", "o");
ylabel("$\Delta f/(B/M)$", "Interpreter", "latex");
yyaxis right;
plot(out, "Marker", "o");
ylim([0, 1]);
grid();
xlim([1, N_SYMBOLS]);
ylabel("SER");
fontsize(10, "points");
xlabel("Symbol");
set(gcf(), "Units", "inches", "Position", [0, 0, 6, 3.75]);
savefig(gcf(), "sat_shift_plot.fig");
exportgraphics(gcf(), "sat_shift_plot.pdf", "ContentType", "vector");