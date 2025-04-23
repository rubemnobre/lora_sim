clear;clc;
addpath("base_functions\");addpath("decider_algorithms\");

n_elevs = 3;
n_snrs = 17;

elevs = linspace(0, pi/2, n_elevs);
snrs = linspace(-14, -6, n_snrs);
freqs = [433e6, 915e6];

fill_table_param(@basic_decider, 10, false, snrs, elevs, freqs, "basic_10.dat");
fill_table_param(@shift_observer2_decider, 10, false, snrs, elevs, freqs, "observer_10.dat");

fill_table_param(@basic_decider, 12, true, snrs, elevs, freqs, "basic_12ldro.dat");
fill_table_param(@shift_observer2_decider, 12, true, snrs, elevs, freqs, "observer_12ldro.dat");

n_elevs = 17;
n_snrs = 3;

elevs = linspace(0, pi/2, n_elevs);
snrs = linspace(-14, -6, n_snrs);
freqs = [433e6, 915e6];

fill_table_param(@basic_decider, 10, false, snrs, elevs, freqs, "basic_10.dat");
fill_table_param(@shift_observer2_decider, 10, false, snrs, elevs, freqs, "observer_10.dat");

fill_table_param(@basic_decider, 12, true, snrs, elevs, freqs, "basic_12ldro.dat");
fill_table_param(@shift_observer2_decider, 12, true, snrs, elevs, freqs, "observer_12ldro.dat");