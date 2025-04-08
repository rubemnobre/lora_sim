function [sequence, symbols] = generate_symbol_sequence(n_symbols, SF, B, OSR, LDRO)
    data = zeros(2^SF*OSR, n_symbols);
    symbols = zeros(1, n_symbols);
    for i = 1:n_symbols
        if LDRO == true
            s = randi([0, 2^(SF-2)-1]);
        else
            s = randi([0, 2^SF-1]);
        end
        symbols(i) = s;
        data(:, i) = generate_lora_symbol(s, SF, B, OSR, LDRO, 0);
    end
    sequence = data(:).'; % apparently transposing conjugates if not using the dot.
end