basic_12ldro = sortrows(readtable("basic_12ldro_sat_link_perf.dat"), "snrs");
basic_10 = sortrows(readtable("basic_10_sat_link_perf.dat"), "snrs");

freq = 433e6;

n_elevs = 5;
elevs = linspace(0, pi/2, n_elevs);
elevs = elevs([1,3,5]);
clf();
hold on;

for i=1:length(elevs)
    basic_12ldro_filtered = basic_12ldro(isapprox(basic_12ldro.elevs, elevs(i), "tight") & basic_12ldro.freqs == freq, :);
    basic_10_filtered = basic_10(isapprox(basic_10.elevs, elevs(i), "tight") & basic_10.freqs == freq, :);

    basic_12ldro_sers = basic_12ldro_filtered.errors./basic_12ldro_filtered.symbols;
    basic_10_sers = basic_10_filtered.errors./basic_10_filtered.symbols;

    s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);

    plot(basic_12ldro_filtered.snrs, basic_12ldro_sers, "DisplayName", "$SF=12$ LDRO; " + s, "LineStyle", "--");
    plot(basic_10_filtered.snrs, basic_10_sers, "DisplayName", "$SF=10$; " + s, "LineStyle", "-");
end

load("basic_curves.mat");

semilogy(snrs_12ldro, sers_12ldro, "DisplayName", "$SF = 12$ LDRO; AWGN only", "LineStyle", "--");
semilogy(snrs_10, sers_10, "DisplayName", "$SF = 10$; AWGN only", "LineStyle", "-");

hold off;
grid();
ylabel("SER");
xlabel("SNR (dB)");
lim = 1e-3;
ylim([lim, 1]);
set(gca, 'YScale', 'log');
fontsize(10, "points");
hl = legend("Location","best", 'FontSize', 9);
set(hl ,'Interpreter','latex');
set(gcf(), "Units", "inches", "Position", [0, 0, 6.25, 3.75]);
ax = gca();
ax.Children = ax.Children([4 3 7 5 1 6 8 2]);
savefig(gcf(), "sat_link_perf_433.fig");
exportgraphics(gcf(), "sat_link_perf_433.pdf", "ContentType", "vector");

freq = 915e6;

n_elevs = 5;
elevs = linspace(0, pi/2, n_elevs);
elevs = elevs([1,3,5]);
clf();
hold on;

for i=1:length(elevs)
    basic_12ldro_filtered = basic_12ldro(isapprox(basic_12ldro.elevs, elevs(i), "tight") & basic_12ldro.freqs == freq, :);
    basic_10_filtered = basic_10(isapprox(basic_10.elevs, elevs(i), "tight") & basic_10.freqs == freq, :);

    basic_12ldro_sers = basic_12ldro_filtered.errors./basic_12ldro_filtered.symbols;
    basic_10_sers = basic_10_filtered.errors./basic_10_filtered.symbols;

    s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);

    plot(basic_12ldro_filtered.snrs, basic_12ldro_sers, "DisplayName", "$SF=12$ LDRO; " + s, "LineStyle", "--");
    plot(basic_10_filtered.snrs, basic_10_sers, "DisplayName", "$SF=10$; " + s, "LineStyle", "-");
end

load("basic_curves.mat");

semilogy(snrs_12ldro, sers_12ldro, "DisplayName", "$SF = 12$ LDRO; AWGN only", "LineStyle", "--");
semilogy(snrs_10, sers_10, "DisplayName", "$SF = 10$; AWGN only", "LineStyle", "-");

hold off;
grid();
ylabel("SER");
xlabel("SNR (dB)");
lim = 1e-3;
ylim([lim, 1]);
set(gca, 'YScale', 'log');
fontsize(10, "points");
hl = legend("Location","best", 'FontSize', 9);
set(hl ,'Interpreter','latex');
set(gcf(), "Units", "inches", "Position", [0, 0, 6.25, 3.75]);
ax = gca();
ax.Children = ax.Children([4 6 3 5 7 1 8 2]);
savefig(gcf(), "sat_link_perf_915.fig");
exportgraphics(gcf(), "sat_link_perf_915.pdf", "ContentType", "vector");