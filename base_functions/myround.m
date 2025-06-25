function y = myround(x)
    y = round(x);
    if abs(x - y) == 0.5
        y = x + randsample([-0.5 0.5], 1);
    end
end

