clear;
clf;
SF = 7;
B = 125000;
GB = 128000;
OSR = 10;
SNR = 100000;
SYMBOL_SHIFT = 32.00;
SHIFT = SYMBOL_SHIFT*B/2^SF;
ZP = 3;

s = 64;

dc = downchirp(SF, B, 1);
Hd = filter_design(B, OSR, GB);
filter_delay = floor(length(Hd.Numerator)/2);
scomp = round(SHIFT*2^SF/B);
scomp_frac = SHIFT*2^SF/B - round(SHIFT*2^SF/B);

ratios = zeros(1, 2^SF);
ratios_flt = zeros(1, 2^SF);
ratios_fltshift = zeros(1, 2^SF);

mag = subplot(2,1,1);
phs = subplot(2,1,2);
hold(mag, "on");
hold(phs, "on");

symbol = generate_lora_symbol(s, SF, B, OSR, 0);
symbol_ds = downsample(symbol, OSR);
dechirped = symbol_ds.*dc;
dechirped_zp = [dechirped zeros(1, (ZP - 1)*length(dechirped))];
sorted = sort(abs(fft(dechirped_zp)), 'descend');
fprintf("largest %d; second %d; ratio %f\n", sorted(1), sorted(2), sorted(2)/sorted(1));
ratios(s) = sorted(2)/sorted(1);
plot(mag, abs(fft(dechirped_zp)), "DisplayName", "No Filter");
plot(phs, unwrap(angle(fft(dechirped_zp))), "DisplayName", "No Filter");

symbol = generate_lora_symbol(s, SF, B, OSR, 0);
symbol_flt_pad = filter(Hd, [symbol zeros(1, filter_delay)]);
symbol_flt = symbol_flt_pad(filter_delay+1:end);
symbol_ds = downsample(symbol_flt, OSR);
dechirped = symbol_ds.*dc;
dechirped_zp = [dechirped zeros(1, (ZP - 1)*length(dechirped))];
spec = fft(dechirped_zp);
sorted = sort(abs(spec), 'descend');
fprintf("largest %d; second %d; ratio %f\n", sorted(1), sorted(2), sorted(2)/sorted(1));
ratios_flt(s) = sorted(2)/sorted(1);
plot(mag, abs(spec), "DisplayName", "Filter");
plot(phs, unwrap(angle(spec)), "DisplayName", "Filter"); 

symbol = generate_lora_symbol(s, SF, B, OSR, SHIFT);
symbol_flt_pad = filter(Hd, [symbol zeros(1, filter_delay)]);
symbol_flt = symbol_flt_pad(filter_delay+1:end);
symbol_ds = downsample(symbol_flt, OSR);
dechirped = symbol_ds.*dc;
dechirped_zp = [dechirped zeros(1, (ZP - 1)*length(dechirped))];
spec = fft(dechirped_zp);
sorted = sort(abs(spec), 'descend');
fprintf("largest %d; second %d; ratio %f\n", sorted(1), sorted(2), sorted(2)/sorted(1));
ratios_fltshift(s) = sorted(2)/sorted(1);
plot(mag, abs(circshift(spec, -scomp*ZP)), "DisplayName", "Shift and Filter");
plot(phs, unwrap(angle(spec)), "DisplayName", "Shift and Filter"); 

linkaxes([mag, phs], "x");
legend();
text = sprintf("Decision FFT Output; Shift = %.2f*B/%d", SYMBOL_SHIFT, 2^SF);
title(text);
xlabel("Bin");
ylabel("Value");
grid();
saveas(gcf(), "compare_plot_fft.png", "png");
hold off;