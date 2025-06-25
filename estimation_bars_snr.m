addpath("base_functions\");addpath("decider_algorithms\");

clf();
hold on;

ZP = 4;
SF = 10;
LDRO = false;
B = 125000;
GB = 128000;
OSR = 2;
SNR = -15;
N = 1000;
N_SHIFTS = 31;
if LDRO == true
    M = 2^SF/4;
else
    M = 2^SF;
end
df = B/M;

shifts = linspace(-0.499*df, 0.499*df, N_SHIFTS);

Hd_filter = filter_design(B, OSR, GB);
Hd = Hd_filter.Numerator;
filter_delay = floor(length(Hd)/2);
dc = downchirp(SF, B, 1);

means = zeros(1, N_SHIFTS);
stdevs = zeros(1, N_SHIFTS);
ests = zeros(1, N);

for j=1:N_SHIFTS
    shift = shifts(j);
    parfor i=1:N
        rng(i);
        s = randi([0, M-1]);
        noise = (sqrt(OSR)/sqrt(2))*randn(1, 2^SF*OSR)*10^(-SNR/20) + (sqrt(OSR)*1j/sqrt(2))*randn(1, 2^SF*OSR)*10^(-SNR/20);
        symbol = generate_lora_symbol(s, SF, B, OSR, LDRO, shift) + noise;
        symbol_flt_pad = filter(Hd, 1, [symbol zeros(1, filter_delay)]);
        symbol_flt = symbol_flt_pad(filter_delay+1:end);
        symbol_ds = downsample(symbol_flt, OSR);
        dechirped = symbol_ds.*dc;
        
        fftres = abs(fft([dechirped zeros(1, length(dechirped)*(ZP-1))]));
        [m, mind] = max(fftres);
        p = (mind-1)/ZP;
        if LDRO == false
            ests(i) = p - myround(p);
        else
            ests(i) = (p/4 - myround(p/4));
        end
    end
    means(j) = mean(ests);
    stdevs(j) = std(ests);
end
st = sprintf("SNR$=%d$; $Z=%d$", SNR, ZP);
errorbar(shifts/df - round(shifts/df), means, stdevs*2, "DisplayName", st);

ZP = 4;
SF = 10;
LDRO = false;
SNR = -19;

Hd_filter = filter_design(B, OSR, GB);
Hd = Hd_filter.Numerator;
filter_delay = floor(length(Hd)/2);
dc = downchirp(SF, B, 1);

means = zeros(1, N_SHIFTS);
stdevs = zeros(1, N_SHIFTS);
ests = zeros(1, N);

for j=1:N_SHIFTS
    shift = shifts(j);
    parfor i=1:N
        rng(i);
        s = randi([0, M-1]);
        noise = (sqrt(OSR)/sqrt(2))*randn(1, 2^SF*OSR)*10^(-SNR/20) + (sqrt(OSR)*1j/sqrt(2))*randn(1, 2^SF*OSR)*10^(-SNR/20);
        symbol = generate_lora_symbol(s, SF, B, OSR, LDRO, shift) + noise;
        symbol_flt_pad = filter(Hd, 1, [symbol zeros(1, filter_delay)]);
        symbol_flt = symbol_flt_pad(filter_delay+1:end);
        symbol_ds = downsample(symbol_flt, OSR);
        dechirped = symbol_ds.*dc;
        
        fftres = abs(fft([dechirped zeros(1, length(dechirped)*(ZP-1))]));
        [m, mind] = max(fftres);
        p = (mind-1)/ZP;
        if LDRO == false
            ests(i) = p - myround(p);
        else
            ests(i) = (p/4 - myround(p/4));
        end
    end
    means(j) = mean(ests);
    stdevs(j) = std(ests);
end

st = sprintf("SNR$=%d$; $Z=%d$", SNR, ZP);
errorbar(shifts/df - round(shifts/df), means, stdevs*2, "DisplayName", st);

% plot(shifts/df - round(shifts/df), shifts/df - round(shifts/df), "DisplayName", "$\Delta f$");
hold off;
grid();
ylabel("$\Delta \hat{\mu}_n$", "Interpreter", "latex");
xlabel("$\Delta f_n$", "Interpreter", "latex");

fontsize(10, "points");
hl = legend("Location", "best", 'FontSize', 9, "Interpreter", "latex");
set(gcf(), "Units", "inches", "Position", [0, 0, 6.25, 3.75 ]);
savefig(gcf(), "estimation_bars_snr.fig");
exportgraphics(gcf(), "estimation_bars_snr.pdf", "ContentType", "vector");