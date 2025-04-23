function [IOI, score] = pickle(tcg, ccurves)
    % shorthand to find best gradient match
    gcge = gradient(ccurves);
    absDiff = abs(tcg - gcge);
    [~, IOI] = min(sum(absDiff, 2));
    score = absDiff(IOI);
end
