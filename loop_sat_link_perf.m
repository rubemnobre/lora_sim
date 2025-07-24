snrs = linspace(-25, -20, 6);

freq = 915e6;

n_elevs = 3;
elevs = linspace(0, 90*pi/180, n_elevs);
lcs = colororder();
clf();
hold on;

% out = sortrows(readtable("ramp_12ldro_sat_link_perf_preamble_5.dat"), "snrs");
% for i=1:length(elevs)
%     out_filtered = out(isapprox(out.elevs, elevs(i), "tight") & out.freqs == freq, :);
%     out_sers = out_filtered.errors./out_filtered.symbols;
%     s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);
%     plot(out_filtered.snrs, out_sers, "DisplayName", "$\rm{SF}=12$ LDRO; " + s + "; Ramp ($\rm{UR}=5$)", "Color", lcs(i,:), "LineStyle", "-.");
% end

out = sortrows(readtable("ramp_12ldro_sat_link_perf_preamble_10.dat"), "snrs");
for i=1:length(elevs)
    out_filtered = out(isapprox(out.elevs, elevs(i), "tight") & out.freqs == freq, :);
    out_sers = out_filtered.errors./out_filtered.symbols;
    s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);
    plot(out_filtered.snrs, out_sers, "DisplayName", "$\rm{SF=12}$ LDRO; " + s + "; Ramp ($\rm{UR}=10$)", "Color", lcs(i,:), "LineStyle","-.");
end

% out = sortrows(readtable("ramp_12ldro_sat_link_perf_preamble_20.dat"), "snrs");
% for i=1:length(elevs)
%     out_filtered = out(isapprox(out.elevs, elevs(i), "tight") & out.freqs == freq, :);
%     out_sers = out_filtered.errors./out_filtered.symbols;
%     s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);
%     plot(out_filtered.snrs, out_sers, "DisplayName", "$\rm{SF=12}$ LDRO; " + s + "; Ramp ($\rm{UR}=2x10$)", "Color", lcs(i,:), "LineStyle","-");
% end

out = sortrows(readtable("observer_12ldro_sat_link_perf_preamble.dat"), "snrs");
for i=1:length(elevs)
    out_filtered = out(isapprox(out.elevs, elevs(i), "tight") & out.freqs == freq, :);
    out_sers = out_filtered.errors./out_filtered.symbols;
    s = sprintf("$E = %.0f^{\\circ}$", elevs(i)*180/pi);
    plot(out_filtered.snrs, out_sers, "DisplayName", "$\rm{SF=12}$ LDRO; " + s + "; Observer", "Color", lcs(i,:), "LineStyle","--");
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