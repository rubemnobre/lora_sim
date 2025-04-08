function res = results(diffs, SF, LDRO)
    if LDRO == true
        M = 2^(SF-2);
    else
        M = 2^(SF);
    end
    diffs = diffs(:);
    res = zeros(1, length(diffs));
    for i=1:length(diffs)
        if diffs(i) > M/2
            res(i) = abs(diffs(i) - M);
        else
            res(i) = diffs(i);
        end
    end
end