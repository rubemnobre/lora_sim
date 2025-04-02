SF = 9;
B = 125000;
GB = 128000;
OSR = 2;
SNR = -9;
N_SYMBOLS = 1000;
Ts = 2^SF/B;
T = N_SYMBOLS*Ts;
ALTITUDE = 500e3;
CENTER_FREQ = 2000e6;
N_SAMPLES = N_SYMBOLS*2^SF*OSR;
ELEVATION = -20*pi/180;
N_RUNS = 20;
N_PAR = 1000;

shifts = satellite_shifts(ALTITUDE, CENTER_FREQ, ELEVATION, B*OSR, N_SAMPLES);
initial_shift = shifts(1);
initial_rate = (shifts(11) - shifts(1))*(10*B*OSR);

Hd_filter = filter_design(B, OSR, GB);
Hd = Hd_filter.Numerator;
filter_delay = floor(length(Hd)/2);

diffs = zeros(N_PAR, N_RUNS, N_SYMBOLS);

[sequence, symbols] = generate_symbol_sequence(N_SYMBOLS, SF, B, OSR);
    
noise = (1/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20) + (1j/sqrt(2))*randn(1, N_SAMPLES)*10^(-SNR/20);
sequence_shifted = shift(sequence, shifts, B*OSR) + noise;
flt_out = filter(Hd, 1, [sequence_shifted zeros(1, filter_delay)]);
sequence_shifted_flt = flt_out(filter_delay+1:end);

r = initial_rate/(B/(Ts*2^SF));

out = shift_observer_decider(sequence_shifted_flt, SF, B, OSR, initial_shift, initial_rate);
sum(abs(sign(out - symbols)))