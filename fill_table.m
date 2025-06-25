clear;clc;
addpath("base_functions\");addpath("decider_algorithms\");

n_elevs = 6;

elevs = linspace(0, 10*pi/180, n_elevs);
freqs = [433e6, 915e6];

% fill_table_param(@basic_decider, 10, false, [-16.5], elevs, freqs, "basic_10_sat_elevs_perf_p.dat");
% 
% fill_table_param(@basic_decider, 12, true, [-21], elevs, freqs, "basic_12ldro_sat_elevs_perf_p.dat");


fill_table_param(@shift_observer2_decider, 10, false, [-16.5], elevs, freqs, "observer_10_sat_elevs_perf_p.dat");

fill_table_param(@shift_observer2_decider, 12, true, [-21], elevs, freqs, "observer_12ldro_sat_elevs_perf_p.dat");