addpath("base_functions\");addpath("decider_algorithms\");
% rng(1);
SF = 12;
LDRO = false;
B = 125000;
GB = 128000;
OSR = 2;
SNR = -22;
N_SYMBOLS = 50;
Ts = 2^SF/B;
T = N_SYMBOLS*Ts;
ALTITUDE = 500e3;
CENTER_FREQ = 915e6;
N_SAMPLES = round((N_SYMBOLS + 12.25)*2^SF*OSR);
ELEVATION = 90*pi/180;
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
initial_rate = (shifts(101) - shifts(1))*(100*B*OSR);
fprintf("rate %f\n", (shifts(2^SF*OSR) - shifts(1))/(B/2^SF));

Hd_filter = filter_design(B, OSR, GB);
Hd = Hd_filter.Numerator;
filter_delay = floor(length(Hd)/2);

diffs = zeros(N_PAR, N_RUNS, N_SYMBOLS);

[sequence, symbols] = generate_symbol_sequence(N_SYMBOLS, SF, B, OSR, LDRO);
preamble = generate_preamble(SF, B, OSR, 8, 52);
sequence = [preamble sequence];
noise = (1/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20) + (1j/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20);
[sequence_shifted, ps] = shiftfun(sequence, shifts, B*OSR);
sequence_shifted = sequence_shifted + noise;
flt_out = filter(Hd, 1, [sequence_shifted zeros(1, filter_delay)]);
sequence_shifted_flt = flt_out(filter_delay+1:end);

r = initial_rate/(B/(Ts*2^SF));

[out, shift_comps] = ramp_decider(sequence_shifted_flt, SF, B, OSR, LDRO, initial_shift, true);


clf();
hold on;
shifts_ds = downsample(shifts(round(12.25*2^SF*OSR)+1:end), 2^SF*OSR)/df;
shifts_out = shift_comps;
plot(shifts_out, "DisplayName", "Estimated");
plot(shifts_ds, "DisplayName", "Actual");
xlabel("Symbol", "FontSize", 10);
ylabel("$\Delta f_n$", "Interpreter", "latex", "FontSize", 10);
hold off;
legend("Location", "northeast", 'FontSize', 9, "Interpreter", "latex");
grid();
fprintf("ser %f\n", sum(abs(sign(out - symbols)))/N_SYMBOLS);
fontsize(10, "points");
% set(gcf(), "Units", "inches", "Position", [0, 0, 6.25, 4.5]);
savefig(gcf(), "loop_perf_12.fig");
exportgraphics(gcf(), "loop_perf_12.pdf", "ContentType", "vector");

% SF = 10;
% LDRO = false;
% SNR = -17.5;
% Ts = 2^SF/B;
% T = N_SYMBOLS*Ts;
% ALTITUDE = 500e3;
% N_SAMPLES = round((N_SYMBOLS + 12.25)*2^SF*OSR);
% N_RUNS = 20;
% N_PAR = 1000;
% if LDRO == true
%     M = 2^SF/4;
% else
%     M = 2^SF;
% end
% df = B/M;
% 
% shifts = satellite_shifts(ALTITUDE, CENTER_FREQ, ELEVATION, B*OSR, N_SAMPLES);
% 
% initial_shift = shifts(1);
% initial_rate = (shifts(11) - shifts(1))*(10*B*OSR);
% 
% Hd_filter = filter_design(B, OSR, GB);
% Hd = Hd_filter.Numerator;
% filter_delay = floor(length(Hd)/2);
% 
% diffs = zeros(N_PAR, N_RUNS, N_SYMBOLS);
% 
% [sequence, symbols] = generate_symbol_sequence(N_SYMBOLS, SF, B, OSR, LDRO);
% preamble = generate_preamble(SF, B, OSR, 8, 52);
% sequence = [preamble sequence];
% noise = (1/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20) + (1j/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20);
% [sequence_shifted, ps] = shiftfun(sequence, shifts, B*OSR);
% sequence_shifted = sequence_shifted + noise;
% flt_out = filter(Hd, 1, [sequence_shifted zeros(1, filter_delay)]);
% sequence_shifted_flt = flt_out(filter_delay+1:end);
% 
% r = initial_rate/(B/(Ts*2^SF));
% 
% [out, shift_comps, res] = shift_observer2_decider(sequence_shifted_flt, SF, B, OSR, LDRO, initial_shift, true);
% 
% clf();
% subplot(211);
% hold on;
% shifts_ds = downsample(shifts(round(12.25*2^SF*OSR)+1:end), 2^SF*OSR)/df;
% shifts_out = shift_comps;
% plot(shifts_out, "DisplayName", "Estimated");
% plot(shifts_ds, "DisplayName", "Actual");
% xlabel("Symbol", "FontSize", 10);
% ylabel("$\Delta f_n$", "Interpreter", "latex", "FontSize", 10);
% hold off;
% legend("Location", "northeast", 'FontSize', 9, "Interpreter", "latex");
% grid();
% subplot(212);
% hold on;
% plot(-res, "DisplayName", "Estimated");
% plot((shifts_out - shifts_ds), "DisplayName", "Actual");
% xlabel("Symbol", "FontSize", 10);
% ylabel("$\mu_n$", "Interpreter", "latex", "FontSize", 10);
% hold off;
% legend("Location", "northeast", 'FontSize', 9, "Interpreter", "latex");
% grid();
% fprintf("ser %f\n", sum(abs(sign(out - symbols)))/N_SYMBOLS);
% fontsize(10, "points");
% set(gcf(), "Units", "inches", "Position", [0, 0, 6.25, 4.5]);
% savefig(gcf(), "loop_perf_10.fig");
% exportgraphics(gcf(), "loop_perf_10.pdf", "ContentType", "vector");