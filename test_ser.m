SF = 7;
B = 125000;
GB = 128000; % Channel center freq increment
OSR = 2;
SNR = -9;
N = 10000;
results = zeros(1, N);

dc = downchirp(SF, B, 1);
Hd_filter = filter_design(B, OSR, GB);
Hd = Hd_filter.Numerator;
filter_delay = floor(length(Hd)/2);

N_SHIFTS = 128;
sers = zeros(1, N_SHIFTS);
bers = zeros(1, N_SHIFTS);
shifts = linspace(0, B/4, N_SHIFTS);

for i=1:N_SHIFTS
    shift = shifts(i);
    scomp = round(shift*2^SF/B);

    parfor j=1:N
        s = randi([0, 2^SF-1]);
        noise = (1/sqrt(2))*randn(1, 2^SF*OSR)*10^(-SNR/20) ...
        + (1j/sqrt(2))*randn(1, 2^SF*OSR)*10^(-SNR/20);

        symbol = generate_lora_symbol(s, SF, B, OSR, shift) + noise;
        symbol_flt_pad = filter(Hd, 1, [symbol zeros(1, filter_delay)]);
        symbol_flt = symbol_flt_pad(filter_delay+1:end);
        symbol_ds = downsample(symbol_flt, OSR);
        dechirped = symbol_ds.*dc;

        fftres = abs(fft(dechirped));

        [m, maxind] = max(fftres);

        err = abs(s - (maxind-1-scomp));

        if err > 2^SF/2
            results(j) = abs(err - 2^SF);
        else
            results(j) = err;
        end
    end
    bers(i) = sum(results)/(N*SF);
    sers(i) = sum(sign(results))/N;
    fprintf("shift %f; ser = %f; ber = %f\n", shifts(i), sers(i), bers(i));
end
subplot(211);
semilogy(shifts*2^SF/B, sers);
ylabel("SER");
xlabel("shift*2^{SF}/B");
ylim([1e-3 1]);
grid();

subplot(212);
semilogy(shifts*2^SF/B, bers);
ylabel("BER");
xlabel("shift*2^{SF}/B");
ylim([1e-3 1]);
grid();

text = sprintf("Error rates vs frequency shift; SF = %d; SNR = %.0f", SF, SNR);
sgtitle(text);

saveas(gcf(), "rates_vs_shift.png", "png");