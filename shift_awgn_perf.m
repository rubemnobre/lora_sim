addpath("base_functions\");
addpath("decider_algorithms\");

clf();
hold on;

SF = 10;
LDRO = false;
B = 125000;
GB = 128000;
OSR = 2;
N = 400000;
results = zeros(1, N);
SHIFT = B/8;

dc = downchirp(SF, B, 1);
Hd_filter = filter_design(B, OSR, GB);
Hd = Hd_filter.Numerator;
filter_delay = floor(length(Hd)/2);

Ts = 2^SF/B;
t = linspace(0, Ts - Ts/2^SF, 2^SF);
scomp = exp(-2j*pi*SHIFT.*t);

N_SNRS = 31;
snrs = linspace(-27, -12, N_SNRS);
sers = zeros(1, N_SNRS);

for i=1:N_SNRS
    snr = snrs(i);

    parfor j=1:N
        rng(j);
        if LDRO == false
            s = randi([0, 2^SF-1]);
        else
            s = randi([0, 2^(SF-2)-1]);
        end
        noise = (sqrt(OSR)/sqrt(2))*randn(1, 2^SF*OSR)*10^(-snr/20) + (sqrt(OSR)*1j/sqrt(2))*randn(1, 2^SF*OSR)*10^(-snr/20);

        symbol = generate_lora_symbol(s, SF, B, OSR, LDRO, SHIFT) + noise;
        symbol_flt_pad = filter(Hd, 1, [symbol zeros(1, filter_delay)]);
        symbol_flt = symbol_flt_pad(filter_delay+1:end);
        symbol_ds = downsample(symbol_flt, OSR);
        dechirped = symbol_ds.*scomp.*dc;

        fftres = abs(fft(dechirped));

        [m, maxind] = max(fftres);
        % fprintf("s=%d; maxind/4=%f\n", s, (maxind-1)/4);

        if LDRO == false
            results(j) = abs(s - (maxind-1));
        else
            results(j) = abs(s - (round((maxind-1)/4)));
        end
    end
    sers(i) = sum(abs(sign(results)))/N;
    fprintf("snr = %f; ser = %f\n", snr, sers(i));
end

st = sprintf("$SF = 10$; $\\Delta f = B/%d$", B/SHIFT);
semilogy(snrs, sers, "DisplayName", st, "LineStyle", "-");

SF = 10;
LDRO = false;
B = 125000;
GB = 128000;
OSR = 2;
results = zeros(1, N);
SHIFT = B/4;

dc = downchirp(SF, B, 1);
Hd_filter = filter_design(B, OSR, GB);
Hd = Hd_filter.Numerator;
filter_delay = floor(length(Hd)/2);

Ts = 2^SF/B;
t = linspace(0, Ts - Ts/2^SF, 2^SF);
scomp = exp(-2j*pi*SHIFT.*t);

N_SNRS = 31;
snrs = linspace(-27, -12, N_SNRS);
sers = zeros(1, N_SNRS);

for i=1:N_SNRS
    snr = snrs(i);

    parfor j=1:N
        rng(j);
        if LDRO == false
            s = randi([0, 2^SF-1]);
        else
            s = randi([0, 2^(SF-2)-1]);
        end
        noise = (sqrt(OSR)/sqrt(2))*randn(1, 2^SF*OSR)*10^(-snr/20) + (sqrt(OSR)*1j/sqrt(2))*randn(1, 2^SF*OSR)*10^(-snr/20);

        symbol = generate_lora_symbol(s, SF, B, OSR, LDRO, SHIFT) + noise;
        symbol_flt_pad = filter(Hd, 1, [symbol zeros(1, filter_delay)]);
        symbol_flt = symbol_flt_pad(filter_delay+1:end);
        symbol_ds = downsample(symbol_flt, OSR);
        dechirped = symbol_ds.*scomp.*dc;

        fftres = abs(fft(dechirped));

        [m, maxind] = max(fftres);
        % fprintf("s=%d; maxind/4=%f\n", s, (maxind-1)/4);

        if LDRO == false
            results(j) = abs(s - (maxind-1));
        else
            results(j) = abs(s - (round((maxind-1)/4)));
        end
    end
    sers(i) = sum(abs(sign(results)))/N;
    fprintf("snr = %f; ser = %f\n", snr, sers(i));
end

st = sprintf("$SF = 10$; $\\Delta f = B/%d$", B/SHIFT);
semilogy(snrs, sers, "DisplayName", st, "LineStyle", "-");

SF = 12;
LDRO = true;
B = 125000;
GB = 128000;
OSR = 2;
results = zeros(1, N);
SHIFT = B/8;

dc = downchirp(SF, B, 1);
Hd_filter = filter_design(B, OSR, GB);
Hd = Hd_filter.Numerator;
filter_delay = floor(length(Hd)/2);

Ts = 2^SF/B;
t = linspace(0, Ts - Ts/2^SF, 2^SF);
scomp = exp(-2j*pi*SHIFT.*t);

N_SNRS = 23;
snrs = linspace(-30, -19, N_SNRS);
sers = zeros(1, N_SNRS);

for i=1:N_SNRS
    snr = snrs(i);

    parfor j=1:N
        rng(j);
        if LDRO == false
            s = randi([0, 2^SF-1]);
        else
            s = randi([0, 2^(SF-2)-1]);
        end
        noise = (sqrt(OSR)/sqrt(2))*randn(1, 2^SF*OSR)*10^(-snr/20) + (sqrt(OSR)*1j/sqrt(2))*randn(1, 2^SF*OSR)*10^(-snr/20);

        symbol = generate_lora_symbol(s, SF, B, OSR, LDRO, SHIFT) + noise;
        symbol_flt_pad = filter(Hd, 1, [symbol zeros(1, filter_delay)]);
        symbol_flt = symbol_flt_pad(filter_delay+1:end);
        symbol_ds = downsample(symbol_flt, OSR);
        dechirped = symbol_ds.*scomp.*dc;

        fftres = abs(fft(dechirped));

        [m, maxind] = max(fftres);
        % fprintf("s=%d; maxind/4=%f\n", s, (maxind-1)/4);

        if LDRO == false
            results(j) = abs(s - (maxind-1));
        else
            results(j) = abs(s - (round((maxind-1)/4)));
        end
    end
    sers(i) = sum(abs(sign(results)))/N;
    fprintf("snr = %f; ser = %f\n", snr, sers(i));
end

st = sprintf("$SF = 12$ LDRO; $\\Delta f = B/%d$", B/SHIFT);
semilogy(snrs, sers, "DisplayName", st, "LineStyle", "--");


SF = 12;
LDRO = true;
B = 125000;
GB = 128000;
OSR = 2;
results = zeros(1, N);
SHIFT = B/4;

dc = downchirp(SF, B, 1);
Hd_filter = filter_design(B, OSR, GB);
Hd = Hd_filter.Numerator;
filter_delay = floor(length(Hd)/2);

Ts = 2^SF/B;
t = linspace(0, Ts - Ts/2^SF, 2^SF);
scomp = exp(-2j*pi*SHIFT.*t);

N_SNRS = 23;
snrs = linspace(-30, -19, N_SNRS);
sers = zeros(1, N_SNRS);

for i=1:N_SNRS
    snr = snrs(i);

    parfor j=1:N
        rng(j);
        if LDRO == false
            s = randi([0, 2^SF-1]);
        else
            s = randi([0, 2^(SF-2)-1]);
        end
        noise = (sqrt(OSR)/sqrt(2))*randn(1, 2^SF*OSR)*10^(-snr/20) + (sqrt(OSR)*1j/sqrt(2))*randn(1, 2^SF*OSR)*10^(-snr/20);

        symbol = generate_lora_symbol(s, SF, B, OSR, LDRO, SHIFT) + noise;
        symbol_flt_pad = filter(Hd, 1, [symbol zeros(1, filter_delay)]);
        symbol_flt = symbol_flt_pad(filter_delay+1:end);
        symbol_ds = downsample(symbol_flt, OSR);
        dechirped = symbol_ds.*scomp.*dc;

        fftres = abs(fft(dechirped));

        [m, maxind] = max(fftres);
        % fprintf("s=%d; maxind/4=%f\n", s, (maxind-1)/4);

        if LDRO == false
            results(j) = abs(s - (maxind-1));
        else
            results(j) = abs(s - (round((maxind-1)/4)));
        end
    end
    sers(i) = sum(abs(sign(results)))/N;
    fprintf("snr = %f; ser = %f\n", snr, sers(i));
end

st = sprintf("$SF = 12$ LDRO; $\\Delta f = B/%d$", B/SHIFT);
semilogy(snrs, sers, "DisplayName", st, "LineStyle", "--");

load("basic_curves.mat");

semilogy(snrs_10, sers_10, "DisplayName", "$SF = 10$; AWGN only");
semilogy(snrs_12ldro, sers_12ldro, "DisplayName", "$SF = 12$ LDRO; AWGN only", "LineStyle", "--");
ylabel("SER");
xlabel("SNR");
ylim([1e-3 1]);
xlim([-30 -13]);
set(gca, 'YScale', 'log');
grid();
hold off;
fontsize(10, "points");
hl = legend("Location","best", 'FontSize', 9);
set(hl ,'Interpreter','latex');
set(gcf(), "Units", "inches", "Position", [0, 0, 6, 4]);
ax = gca();
ax.Children = ax.Children([5 6 2 3 4 1]);
savefig(gcf(), "shift_awgn_perf.fig");
exportgraphics(gcf(), "shift_awgn_perf.pdf", "ContentType", "vector");


fontsize(10, "points");
hl = legend("Location","best", 'FontSize', 7);
set(hl ,'Interpreter','latex');
set(gcf(), "Units", "inches", "Position", [0, 0, 4.2, 3]);
exportgraphics(gcf(), "shift_awgn_perf_pres.pdf", "ContentType", "vector");