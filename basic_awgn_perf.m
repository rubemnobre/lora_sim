addpath("base_functions\");
addpath("decider_algorithms\");

clf();
hold on;

SF = 12;
LDRO = true;
B = 125000;
OSR = 1;
N = 400000;
results = zeros(1, N);

dc = downchirp(SF, B, 1);

N_SNRS = 19;
snrs = linspace(-30, -21, N_SNRS);
sers = zeros(1, N_SNRS);

for i=1:N_SNRS
    shift = 0;
    scomp = round(shift*2^SF/B);
    snr = snrs(i);

    parfor j=1:N
        rng(j);
        if LDRO == false
            s = randi([0, 2^SF-1]);
        else
            s = randi([0, 2^(SF-2)-1]);
        end
        noise = (sqrt(OSR)/sqrt(2))*randn(1, 2^SF*OSR)*10^(-snr/20) + (sqrt(OSR)*1j/sqrt(2))*randn(1, 2^SF*OSR)*10^(-snr/20);

        symbol = generate_lora_symbol(s, SF, B, OSR, LDRO, shift) + noise;
        symbol_ds = downsample(symbol, OSR);
        dechirped = symbol_ds.*dc;

        fftres = abs(fft(dechirped));

        [m, maxind] = max(fftres);
        % fprintf("s=%d; maxind/4=%f\n", s, (maxind-1)/4);

        if LDRO == false
            results(j) = abs(s - (maxind-1-scomp));
        else
            results(j) = abs(s - (round((maxind-1)/4)-scomp));
        end
    end
    sers(i) = sum(abs(sign(results)))/N;
    fprintf("snr = %f; ser = %f\n", snr, sers(i));
end
snrs_12ldro = snrs;
sers_12ldro = sers;

semilogy(snrs, sers, "DisplayName", "$SF = 12$ with LDRO", "LineStyle","--");

SF = 10;
LDRO = false;
B = 125000;
OSR = 1;
N = 400000;
results = zeros(1, N);

dc = downchirp(SF, B, 1);

N_SNRS = 23;
snrs = linspace(-27, -16, N_SNRS);
sers = zeros(1, N_SNRS);

for i=1:N_SNRS
    shift = 0;
    scomp = round(shift*2^SF/B);
    snr = snrs(i);

    parfor j=1:N
        rng(j);
        if LDRO == false
            s = randi([0, 2^SF-1]);
        else
            s = randi([0, 2^(SF-2)-1]);
        end
        noise = (sqrt(OSR)/sqrt(2))*randn(1, 2^SF*OSR)*10^(-snr/20) + (sqrt(OSR)*1j/sqrt(2))*randn(1, 2^SF*OSR)*10^(-snr/20);

        symbol = generate_lora_symbol(s, SF, B, OSR, LDRO, shift) + noise;
        symbol_ds = downsample(symbol, OSR);
        dechirped = symbol_ds.*dc;

        fftres = abs(fft(dechirped));

        [m, maxind] = max(fftres);
        % fprintf("s=%d; maxind/4=%f\n", s, (maxind-1)/4);

        if LDRO == false
            results(j) = abs(s - (maxind-1-scomp));
        else
            results(j) = abs(s - (round((maxind-1)/4)-scomp));
        end
    end
    sers(i) = sum(abs(sign(results)))/N;
    fprintf("snr = %f; ser = %f\n", snr, sers(i));
end
snrs_10 = snrs;
sers_10 = sers;
semilogy(snrs, sers, "DisplayName", "$SF = 10$", "LineStyle", "-");

ylabel("SER");
xlabel("SNR (dB)");
ylim([1e-3 1]);
set(gca, 'YScale', 'log');
grid();
hold off;
fontsize(10, "points");
hl = legend("Location","best", 'FontSize', 9);
set(hl ,'Interpreter','latex');
set(gcf(), "Units", "inches", "Position", [0, 0, 4.5, 3]);
savefig(gcf(), "basic_awgn_perf.fig");
exportgraphics(gcf(), "basic_awgn_perf.pdf", "ContentType", "vector");
save("basic_curves.mat","sers_10", "snrs_10", "sers_12ldro", "snrs_12ldro");