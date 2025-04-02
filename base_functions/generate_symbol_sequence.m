function [sequence, symbols] = generate_symbol_sequence(n_symbols, SF, B, OSR)
    data = zeros(2^SF*OSR, n_symbols);
    symbols = zeros(1, n_symbols);
    for i = 1:n_symbols
        s  = randi([0, 2^SF-1]);
        symbols(i) = s;
        data(:, i) = generate_lora_symbol(s, SF, B, OSR, 0);
    end
    sequence = data(:).'; % apparently transposing conjugates without the dot.
end