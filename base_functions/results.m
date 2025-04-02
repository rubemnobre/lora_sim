function res = results(diffs, SF)
    diffs = diffs(:);
    res = zeros(1, length(diffs));
    for i=1:length(diffs)
        if diffs(i) > 2^SF/2
            res(i) = abs(diffs(i) - 2^SF);
        else
            res(i) = diffs(i);
        end
    end
end