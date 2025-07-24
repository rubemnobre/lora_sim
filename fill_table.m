clear;clc;
addpath("base_functions\");addpath("decider_algorithms\");

n_elevs = 3;
snrs = linspace(-25, -20, 6);
elevs = linspace(0, 90*pi/180, n_elevs);
freqs = [915e6];

% fill_table_param(@ramp_decider, 12, true, snrs, elevs, freqs, "ramp_12ldro_sat_link_perf_preamble_10.dat");

fill_table_param(@shift_observer2_decider, 12, true, snrs, elevs, freqs, "observer_12ldro_sat_link_perf_preamble.dat");