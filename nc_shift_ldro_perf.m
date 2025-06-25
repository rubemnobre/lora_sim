addpath("base_functions\");
addpath("decider_algorithms\");

clf();
hold on;

SF = 12;
LDRO = false;
B = 125000;
GB = 128000;
OSR = 1;
N = 200000;
results = zeros(1, N);
SHIFT = (B/2^SF)/2;

dc = downchirp(SF, B, 1);

Ts = 2^SF/B;
t = linspace(0, Ts - Ts/2^SF, 2^SF);

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
        dechirped = symbol.*dc;

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

st = sprintf("SF $= 12$; $\\Delta f = (B/M)/%d$", (B/2^SF)/SHIFT);
semilogy(snrs, sers, "DisplayName", st, "LineStyle", "--");

SF = 12;
LDRO = false;
B = 125000;
GB = 128000;
OSR = 1;
results = zeros(1, N);
SHIFT = (B/2^SF)/4;

dc = downchirp(SF, B, 1);

Ts = 2^SF/B;
t = linspace(0, Ts - Ts/2^SF, 2^SF);

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
        dechirped = symbol.*dc;

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

st = sprintf("SF $= 12$; $\\Delta f = (B/2^{SF})/%d$", (B/2^SF)/SHIFT);
semilogy(snrs, sers, "DisplayName", st, "LineStyle", "--");

SF = 12;
LDRO = true;
B = 125000;
GB = 128000;
OSR = 1;
results = zeros(1, N);
SHIFT = (B/2^SF)/2;

dc = downchirp(SF, B, 1);

Ts = 2^SF/B;
t = linspace(0, Ts - Ts/2^SF, 2^SF);

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
        dechirped = symbol.*dc;

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

st = sprintf("SF $= 12$ LDRO; $\\Delta f = (B/2^{SF})/%d$", (B/2^SF)/SHIFT);
semilogy(snrs, sers, "DisplayName", st, "LineStyle", "--");


SF = 12;
LDRO = true;
B = 125000;
GB = 128000;
OSR = 1;
results = zeros(1, N);
SHIFT = (B/2^SF)/4;

dc = downchirp(SF, B, 1);

Ts = 2^SF/B;
t = linspace(0, Ts - Ts/2^SF, 2^SF);

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
        dechirped = symbol.*dc;

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

st = sprintf("SF $= 12$ LDRO; $\\Delta f = (B/2^{SF})/%d$", (B/2^SF)/SHIFT);
semilogy(snrs, sers, "DisplayName", st, "LineStyle", "--");

load("basic_curves.mat");
semilogy(snrs_12ldro, sers_12ldro, "DisplayName", "$SF = 12$ LDRO");
ylabel("SER");
xlabel("SNR (dB)");
ylim([1e-3 1]);
xlim([-30 -19]);
set(gca, 'YScale', 'log');
grid();
hold off;
ax = gca();
ax.Children = ax.Children([5 4 2 3 1]);

fontsize(10, "points");
hl = legend("Location","best", 'FontSize', 9);
set(hl ,'Interpreter','latex');
set(gcf(), "Units", "inches", "Position", [0, 0, 6, 4]);
savefig(gcf(), "nc_shift_ldro_perf.fig");
exportgraphics(gcf(), "nc_shift_ldro_perf.pdf", "ContentType", "vector");


fontsize(10, "points");
hl = legend("Location","best", 'FontSize', 7);
set(hl ,'Interpreter','latex');
set(gcf(), "Units", "inches", "Position", [0, 0, 3.8, 3*3.8/4.2]);
exportgraphics(gcf(), "nc_shift_ldro_perf_pres.pdf", "ContentType", "vector");