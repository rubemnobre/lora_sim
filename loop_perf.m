addpath("base_functions\");addpath("decider_algorithms\");
% rng(1);
SF = 12;
LDRO = true;
B = 125000;
GB = 128000;
OSR = 2;
SNR = -23;
N_SYMBOLS = 100;
Ts = 2^SF/B;
T = N_SYMBOLS*Ts;
ALTITUDE = 500e3;
CENTER_FREQ = 915e6;
N_SAMPLES = N_SYMBOLS*2^SF*OSR;
ELEVATION = 63*pi/180;
N_RUNS = 20;
N_PAR = 1000;
if LDRO == true
    M = 2^SF/4;
else
    M = 2^SF;
end
df = B/M;

shifts = satellite_shifts(ALTITUDE, CENTER_FREQ, ELEVATION, B*OSR, N_SAMPLES);

initial_shift = shifts(1);
initial_rate = (shifts(11) - shifts(1))*(10*B*OSR);

Hd_filter = filter_design(B, OSR, GB);
Hd = Hd_filter.Numerator;
filter_delay = floor(length(Hd)/2);

diffs = zeros(N_PAR, N_RUNS, N_SYMBOLS);

[sequence, symbols] = generate_symbol_sequence(N_SYMBOLS, SF, B, OSR, LDRO);
    
noise = (1/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20) + (1j/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20);
[sequence_shifted, ps] = shiftfun(sequence, shifts, B*OSR);
sequence_shifted = sequence_shifted + noise;
flt_out = filter(Hd, 1, [sequence_shifted zeros(1, filter_delay)]);
sequence_shifted_flt = flt_out(filter_delay+1:end);

r = initial_rate/(B/(Ts*2^SF));

[out, shift_comps, res] = shift_observer2_decider(sequence_shifted_flt, SF, B, OSR, LDRO, initial_shift, initial_rate);


clf();
subplot(211);
hold on;
shifts_ds = downsample(shifts, 2^SF*OSR)/df;
shifts_out = shift_comps;
plot(shifts_out, "DisplayName", "Estimated");
plot(shifts_ds, "DisplayName", "Actual");
xlabel("Symbol", "FontSize", 10);
ylabel("$\Delta f_n$", "Interpreter", "latex", "FontSize", 10);
hold off;
legend("Location", "northeast", 'FontSize', 9, "Interpreter", "latex");
grid();
subplot(212);
hold on;
plot(-res, "DisplayName", "Estimated");
plot((shifts_out - shifts_ds), "DisplayName", "Actual");
xlabel("Symbol", "FontSize", 10);
ylabel("$\mu_n$", "Interpreter", "latex", "FontSize", 10);
hold off;
legend("Location", "northeast", 'FontSize', 9, "Interpreter", "latex");
grid();
fprintf("ser %f\n", sum(abs(sign(out - symbols)))/N_SYMBOLS);
fontsize(10, "points");
set(gcf(), "Units", "inches", "Position", [0, 0, 6.25, 4.5]);
savefig(gcf(), "loop_perf_12.fig");
exportgraphics(gcf(), "loop_perf_12.pdf", "ContentType", "vector");

SF = 10;
LDRO = false;
SNR = -17.5;
Ts = 2^SF/B;
T = N_SYMBOLS*Ts;
ALTITUDE = 500e3;
N_SAMPLES = N_SYMBOLS*2^SF*OSR;
N_RUNS = 20;
N_PAR = 1000;
if LDRO == true
    M = 2^SF/4;
else
    M = 2^SF;
end
df = B/M;

shifts = satellite_shifts(ALTITUDE, CENTER_FREQ, ELEVATION, B*OSR, N_SAMPLES);

initial_shift = shifts(1);
initial_rate = (shifts(11) - shifts(1))*(10*B*OSR);

Hd_filter = filter_design(B, OSR, GB);
Hd = Hd_filter.Numerator;
filter_delay = floor(length(Hd)/2);

diffs = zeros(N_PAR, N_RUNS, N_SYMBOLS);

[sequence, symbols] = generate_symbol_sequence(N_SYMBOLS, SF, B, OSR, LDRO);
    
noise = (1/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20) + (1j/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20);
[sequence_shifted, ps] = shiftfun(sequence, shifts, B*OSR);
sequence_shifted = sequence_shifted + noise;
flt_out = filter(Hd, 1, [sequence_shifted zeros(1, filter_delay)]);
sequence_shifted_flt = flt_out(filter_delay+1:end);

r = initial_rate/(B/(Ts*2^SF));

[out, shift_comps, res] = shift_observer2_decider(sequence_shifted_flt, SF, B, OSR, LDRO, initial_shift, initial_rate);

clf();
subplot(211);
hold on;
shifts_ds = downsample(shifts, 2^SF*OSR)/df;
shifts_out = shift_comps;
plot(shifts_out, "DisplayName", "Estimated");
plot(shifts_ds, "DisplayName", "Actual");
xlabel("Symbol", "FontSize", 10);
ylabel("$\Delta f_n$", "Interpreter", "latex", "FontSize", 10);
hold off;
legend("Location", "northeast", 'FontSize', 9, "Interpreter", "latex");
grid();
subplot(212);
hold on;
plot(-res, "DisplayName", "Estimated");
plot((shifts_out - shifts_ds), "DisplayName", "Actual");
xlabel("Symbol", "FontSize", 10);
ylabel("$\mu_n$", "Interpreter", "latex", "FontSize", 10);
hold off;
legend("Location", "northeast", 'FontSize', 9, "Interpreter", "latex");
grid();
fprintf("ser %f\n", sum(abs(sign(out - symbols)))/N_SYMBOLS);
fontsize(10, "points");
set(gcf(), "Units", "inches", "Position", [0, 0, 6.25, 4.5]);
savefig(gcf(), "loop_perf_10.fig");
exportgraphics(gcf(), "loop_perf_10.pdf", "ContentType", "vector");