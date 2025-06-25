function fill_table_param(decider, sf, ldro, snrs, elevs, freqs, path)
    
    if isfile(path) == false
        elevs = [];
        snrs = [];
        freqs = [];
        errors = [];
        symbols = [];
        t = table(elevs, snrs, freqs, errors, symbols);
        writetable(t, path);
    end
    
    configs = combinations(elevs, snrs, freqs);

    for j=1:height(configs)
        t = readtable(path);
        elev = configs.elevs(j);
        snr = configs.snrs(j);
        freq = configs.freqs(j);
        fprintf("%s - running config %d/%d. elev = %f; snr = %f, sf = %d, ldro = %d, freq = %f\n", path, j, height(configs), elev, snr, sf, ldro, freq);
        if height(t) == 0
            [err, symbols] = sat_test_param(decider, elev, snr, sf, ldro, freq);
            t = [t; {elev, snr, freq, err, symbols}];
            writetable(t, path);
        elseif isempty(find(round(t.elevs, 3) == round(elev, 3) & t.freqs == freq & t.snrs == snr, 1)) == true
            [err, symbols] = sat_test_param(decider, elev, snr, sf, ldro, freq);
            t = [t; {elev, snr, freq, err, symbols}];
            writetable(t, path);
        else
            fprintf("skipping, point already in table\n");
        end
    end
end