observer_12ldro = sortrows(readtable("observer_12ldro_sat_elevs_perf_p.dat"), "elevs");
basic_12ldro = sortrows(readtable("basic_12ldro_sat_elevs_perf_p.dat"), "elevs");
lcs = colororder();
freq = 433e6;

n_elevs = 19;
elevs = linspace(0, pi/2, n_elevs);
snrs = [-21];

clf();
hold on;

for i=1:length(snrs)
    basic_12ldro_filtered = basic_12ldro(isapprox(basic_12ldro.snrs, snrs(i), "tight") & basic_12ldro.freqs == freq, :);
    basic_12ldro_sers = basic_12ldro_filtered.errors./basic_12ldro_filtered.symbols;

    observer_12ldro_filtered = observer_12ldro(isapprox(observer_12ldro.snrs, snrs(i), "tight") & observer_12ldro.freqs == freq, :);
    observer_12ldro_sers = observer_12ldro_filtered.errors./observer_12ldro_filtered.symbols;

    s = sprintf("$SNR=%.1f$ dB", snrs(i));

    plot(basic_12ldro_filtered.elevs.*180/pi, basic_12ldro_sers, "DisplayName", "$f_c=433$ MHz " + s + " Static", "Color", lcs(1,:), "LineStyle", "--");
    plot(observer_12ldro_filtered.elevs.*180/pi, observer_12ldro_sers, "DisplayName", "$f_c=433$ MHz " + s + " Proposed", "Color", lcs(1,:), "LineStyle", "-.");
end

freq = 915e6;

n_elevs = 19;
elevs = linspace(0, pi/2, n_elevs);
snrs = [-21];

for i=1:length(snrs)
    basic_12ldro_filtered = basic_12ldro(isapprox(basic_12ldro.snrs, snrs(i), "tight") & basic_12ldro.freqs == freq, :);
    basic_12ldro_sers = basic_12ldro_filtered.errors./basic_12ldro_filtered.symbols;

    observer_12ldro_filtered = observer_12ldro(isapprox(observer_12ldro.snrs, snrs(i), "tight") & observer_12ldro.freqs == freq, :);
    observer_12ldro_sers = observer_12ldro_filtered.errors./observer_12ldro_filtered.symbols;

    s = sprintf("$SNR=%.1f$ dB", snrs(i));

    plot(basic_12ldro_filtered.elevs.*180/pi, basic_12ldro_sers, "DisplayName", "$f_c=915$ MHz " + s + " Static", "Color", lcs(2,:), "LineStyle", "--");
    plot(observer_12ldro_filtered.elevs.*180/pi, observer_12ldro_sers, "DisplayName", "$f_c=915$ MHz " + s + " Proposed", "Color", lcs(2,:), "LineStyle", "-.");
end


hold off;
grid();
ylabel("SER");
xlabel("E");
lim = 1e-3;
ylim([lim, 1]);
set(gca, 'YScale', 'log');
fontsize(10, "points");
hl = legend("Location","best", 'FontSize', 9);
set(hl ,'Interpreter','latex');
set(gcf(), "Units", "inches", "Position", [0, 0, 6.25, 3.75]);
ax = gca();

savefig(gcf(), "loop_sat_elevs_perf_12.fig");
exportgraphics(gcf(), "loop_sat_elevs_perf_12.pdf", "ContentType", "vector");



observer_10 = sortrows(readtable("observer_10_sat_elevs_perf_p.dat"), "elevs");
basic_10 = sortrows(readtable("basic_10_sat_elevs_perf_p.dat"), "elevs");

freq = 433e6;

n_elevs = 19;
elevs = linspace(0, pi/2, n_elevs);
snrs = [-16.5];

clf();
hold on;

for i=1:length(snrs)
    basic_10_filtered = basic_10(isapprox(basic_10.snrs, snrs(i), "tight") & basic_10.freqs == freq, :);
    basic_10_sers = basic_10_filtered.errors./basic_10_filtered.symbols;

    observer_10_filtered = observer_10(isapprox(observer_10.snrs, snrs(i), "tight") & observer_10.freqs == freq, :);
    observer_10_sers = observer_10_filtered.errors./observer_10_filtered.symbols;

    s = sprintf("$SNR=%.1f$ dB", snrs(i));

    plot(basic_10_filtered.elevs.*180/pi, basic_10_sers, "DisplayName", "$f_c=433$ MHz " + s + " Static", "Color", lcs(1,:), "LineStyle", "-");
    plot(observer_10_filtered.elevs.*180/pi, observer_10_sers, "DisplayName", "$f_c=433$ MHz " + s + " Proposed", "Color", lcs(1,:), "LineStyle", ":", "LineWidth", 1.5);
end

freq = 915e6;

n_elevs = 19;
elevs = linspace(0, pi/2, n_elevs);
snrs = [-16.5];

for i=1:length(snrs)
    basic_10_filtered = basic_10(isapprox(basic_10.snrs, snrs(i), "tight") & basic_10.freqs == freq, :);
    basic_10_sers = basic_10_filtered.errors./basic_10_filtered.symbols;

    observer_10_filtered = observer_10(isapprox(observer_10.snrs, snrs(i), "tight") & observer_10.freqs == freq, :);
    observer_10_sers = observer_10_filtered.errors./observer_10_filtered.symbols;

    s = sprintf("$SNR=%.1f$ dB", snrs(i));

    plot(basic_10_filtered.elevs.*180/pi, basic_10_sers, "DisplayName", "$f_c=915$ MHz " + s + " Static", "Color", lcs(2,:), "LineStyle", "-");
    plot(observer_10_filtered.elevs.*180/pi, observer_10_sers, "DisplayName", "$f_c=915$ MHz " + s + " Proposed", "Color", lcs(2,:), "LineStyle", ":", "LineWidth", 1.5);
end


hold off;
grid();
ylabel("SER");
xlabel("E");
lim = 1e-3;
ylim([lim, 1]);
xlim([0, 90]);
set(gca, 'YScale', 'log');
fontsize(10, "points");
hl = legend("Location","best", 'FontSize', 9);
set(hl ,'Interpreter','latex');
set(gcf(), "Units", "inches", "Position", [0, 0, 6.25, 3.75]);
ax = gca();

savefig(gcf(), "loop_sat_elevs_perf_10.fig");
exportgraphics(gcf(), "loop_sat_elevs_perf_10.pdf", "ContentType", "vector");