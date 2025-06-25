freq = 433e6;

n_elevs = 5;
elevs = linspace(0, pi/2, n_elevs);
elevs = elevs([1,3,5]);
lcs = colororder();
clf();
hold on;

out = sortrows(readtable("observer_12ldro_sat_link_perf.dat"), "snrs");
for i=1:length(elevs)
    out_filtered = out(isapprox(out.elevs, elevs(i), "tight") & out.freqs == freq, :);
    out_sers = out_filtered.errors./out_filtered.symbols;
    s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);
    plot(out_filtered.snrs, out_sers, "DisplayName", "$SF=12$ LDRO; " + s + "; Proposed", "Color", lcs(i,:), "LineStyle", "-.");
end

out = sortrows(readtable("basic_12ldro_sat_link_perf.dat"), "snrs");
for i=1:length(elevs)
    out_filtered = out(isapprox(out.elevs, elevs(i), "tight") & out.freqs == freq, :);
    out_sers = out_filtered.errors./out_filtered.symbols;
    s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);
    plot(out_filtered.snrs, out_sers, "DisplayName", "$SF=12$ LDRO; " + s + "; Static", "Color", lcs(i,:), "LineStyle","--");
end

out = sortrows(readtable("observer_10_sat_link_perf.dat"), "snrs");
for i=1:length(elevs)
    out_filtered = out(isapprox(out.elevs, elevs(i), "tight") & out.freqs == freq, :);
    out_sers = out_filtered.errors./out_filtered.symbols;
    s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);
    plot(out_filtered.snrs, out_sers, "DisplayName", "$SF=10$; " + s + "; Proposed", "Color", lcs(i,:), "LineStyle", ":", "LineWidth", 1.5);
end

out = sortrows(readtable("basic_10_sat_link_perf.dat"), "snrs");
for i=1:length(elevs)
    out_filtered = out(isapprox(out.elevs, elevs(i), "tight") & out.freqs == freq, :);
    out_sers = out_filtered.errors./out_filtered.symbols;
    s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);
    plot(out_filtered.snrs, out_sers, "DisplayName", "$SF=10$; " + s + "; Static", "Color", lcs(i,:), "LineStyle","-");
end

hold off;
grid();
ylabel("SER");
xlabel("SNR (dB)");
lim = 1e-3;
ylim([lim, 1]);
set(gca, 'YScale', 'log');
fontsize(10, "points");
hl = legend("Location","best", 'FontSize', 8);
set(hl ,'Interpreter','latex');
set(gcf(), "Units", "inches", "Position", [0, 0, 6.25, 3.75]);
ax = gca();

savefig(gcf(), "loop_sat_link_perf_433.fig");
exportgraphics(gcf(), "loop_sat_link_perf_433.pdf", "ContentType", "vector");



freq = 915e6;

n_elevs = 5;
elevs = linspace(0, pi/2, n_elevs);
elevs = elevs([1,3,5]);
lcs = colororder();
clf();
hold on;

out = sortrows(readtable("observer_12ldro_sat_link_perf.dat"), "snrs");
for i=1:length(elevs)
    out_filtered = out(isapprox(out.elevs, elevs(i), "tight") & out.freqs == freq, :);
    out_sers = out_filtered.errors./out_filtered.symbols;
    s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);
    plot(out_filtered.snrs, out_sers, "DisplayName", "$SF=12$ LDRO; " + s + "; Proposed", "Color", lcs(i,:), "LineStyle", "-.");
end

out = sortrows(readtable("basic_12ldro_sat_link_perf.dat"), "snrs");
for i=1:length(elevs)
    out_filtered = out(isapprox(out.elevs, elevs(i), "tight") & out.freqs == freq, :);
    out_sers = out_filtered.errors./out_filtered.symbols;
    s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);
    plot(out_filtered.snrs, out_sers, "DisplayName", "$SF=12$ LDRO; " + s + "; Static", "Color", lcs(i,:), "LineStyle","--");
end

out = sortrows(readtable("observer_10_sat_link_perf.dat"), "snrs");
for i=1:length(elevs)
    out_filtered = out(isapprox(out.elevs, elevs(i), "tight") & out.freqs == freq, :);
    out_sers = out_filtered.errors./out_filtered.symbols;
    s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);
    plot(out_filtered.snrs, out_sers, "DisplayName", "$SF=10$; " + s + "; Proposed", "Color", lcs(i,:), "LineStyle", ":", "LineWidth", 1.5);
end

out = sortrows(readtable("basic_10_sat_link_perf.dat"), "snrs");
for i=1:length(elevs)
    out_filtered = out(isapprox(out.elevs, elevs(i), "tight") & out.freqs == freq, :);
    out_sers = out_filtered.errors./out_filtered.symbols;
    s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);
    plot(out_filtered.snrs, out_sers, "DisplayName", "$SF=10$; " + s + "; Static", "Color", lcs(i,:), "LineStyle","-");
end

hold off;
grid();
ylabel("SER");
xlabel("SNR (dB)");
lim = 1e-3;
ylim([lim, 1]);
set(gca, 'YScale', 'log');
fontsize(10, "points");
hl = legend("Location","best", 'FontSize', 8);
set(hl ,'Interpreter','latex');
set(gcf(), "Units", "inches", "Position", [0, 0, 6.25, 3.75]);
ax = gca();

savefig(gcf(), "loop_sat_link_perf_915.fig");
exportgraphics(gcf(), "loop_sat_link_perf_915.pdf", "ContentType", "vector");