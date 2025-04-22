addpath("base_functions\");addpath("decider_algorithms\");

SF = 12;
LDRO = true;
B = 125000;
GB = 128000;
OSR = 2;
SNR = -24;
N_SYMBOLS = 50;
Ts = 2^SF/B;
T = N_SYMBOLS*Ts;
ALTITUDE = 500e3;
CENTER_FREQ = 915e6;
N_SAMPLES = N_SYMBOLS*2^SF*OSR;
ELEVATION = 0*pi/180;
N_RUNS = 20;
N_PAR = 1000;

shifts = satellite_shifts(ALTITUDE, CENTER_FREQ, ELEVATION, B*OSR, N_SAMPLES);
initial_shift = shifts(1);
initial_rate = (shifts(11) - shifts(1))*(10*B*OSR);

Hd_filter = filter_design(B, OSR, GB);
Hd = Hd_filter.Numerator;
filter_delay = floor(length(Hd)/2);

diffs = zeros(N_PAR, N_RUNS, N_SYMBOLS);

[sequence, symbols] = generate_symbol_sequence(N_SYMBOLS, SF, B, OSR, LDRO);
    
noise = (1/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20) + (1j/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20);
sequence_shifted = shift(sequence, shifts, B*OSR) + noise;
flt_out = filter(Hd, 1, [sequence_shifted zeros(1, filter_delay)]);
sequence_shifted_flt = flt_out(filter_delay+1:end);

r = initial_rate/(B/(Ts*2^SF));

[out, shift_comps, res] = shift_observer2_decider(sequence_shifted_flt, SF, B, OSR, LDRO, initial_shift, initial_rate);
clf();
subplot(211);
hold on;
shifts_ds = downsample(shifts, 2^SF*OSR);
shifts_out = (-shift_comps/2 + -shift_comps(1)/2)*B/2^SF; %% formula is kind of ad-hoc
plot(shifts_ds, "DisplayName", "Doppler Shift");
plot(shifts_out, "DisplayName", "Observer Compensation");
xlabel("Symbol");
ylabel("Frequency (Hz)");
hold off;
legend();
grid();
subplot(212);
hold on;
plot(0.5*res*B/2^SF, "DisplayName", "Estimated Residue");
plot(shifts_out - shifts_ds, "DisplayName", "Actual Residue");
xlabel("Symbol");
ylabel("Frequency (Hz)");
hold off;
legend();
grid();

t = sprintf("Satellite Altitude = %.0f km; Center Frequency = %.0f MHz; \nSF = %d; Elevation = %.1f deg; BW = %.0f kHz; SNR = %d dB", ALTITUDE/1e3, CENTER_FREQ/1e6, SF, ELEVATION*180/pi, B/1e3, SNR);
sgtitle(t);

if LDRO == true
    ldro_str = "ldro";
else
    ldro_str = "";
end
t = sprintf("residue_%d_%d_%d_%d_%s", round(CENTER_FREQ/1e6), ELEVATION*180/pi, SF, SNR, ldro_str);
saveas(gcf(), t + ".png", "png");
save(t + ".mat");