if isfile("observer.dat") == false
    elevs = [];
    snrs = [];
    freqs = [];
    errors = [];
    symbols = [];
    t = table(elevs, snrs, freqs, errors, symbols);
    writetable(t, "observer.dat");
end
if isfile("basic.dat") == false
    elevs = [];
    snrs = [];
    freqs = [];
    errors = [];
    symbols = [];
    t = table(elevs, snrs, freqs, errors, symbols);
    writetable(t, "basic.dat");
end

n_elevs = 21;
n_snrs = 5;
t = readtable("observer.dat");

elevs = linspace(0, pi/2, n_elevs);
snrs = linspace(-12, -6, n_snrs);

for j=1:n_snrs
    for i=1:n_elevs
        elev = elevs(i);
        fprintf("elev %f\n", elev*180/pi);
        sf = 10;
        snr = snrs(j) - 3*(sf-7);
        freq = 433e6;
        if height(t) == 0
            [error, symbols] = sat_test_param(@shift_observer2_decider, elev, snr, sf, false, freq);
        elseif isempty(find(t.elevs == elev & t.freqs == freq & t.snrs == snr)) == true
            [error, symbols] = sat_test_param(@shift_observer2_decider, elev, snr, 10, false, freq);
        end
        t = [t; {elev, snr, freq, error, symbols}];
        writetable(t, "observer.dat");
    end
    
    
    t = readtable("basic.dat");
    
    for i=1:n_elevs
        elev = elevs(i);
        fprintf("elev %f\n", elev*180/pi);
        sf = 10;
        snr = snrs(j) - 3*(sf-7);
        freq = 433e6;
        if height(t) == 0
            [error, symbols] = sat_test_param(@basic_decider, elev, snr, sf, false, freq);
        elseif isempty(find(t.elevs == elev & t.freqs == freq & t.snrs == snr)) == true
            [error, symbols] = sat_test_param(@basic_decider, elev, snr, 10, false, freq);
        end
        t = [t; {elev, snr, freq, error, symbols}];
        writetable(t, "basic.dat");
    end
end