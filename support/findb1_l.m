function b1 = findb1_l(R, Aa, Bb, lim)
    b1 = zeros(size(R)); 
    % Case 1: R <= lim
    idx1 = R <= lim;
    b1(idx1) = Aa .* sqrt(Bb - (R(idx1) ./ 2).^2);

    % Case 2: lim < R <= 2*lim
    idx2 = R > lim & R <= 2 * lim;
    b1(idx2) = Aa .* sqrt((4 * Bb - (R(idx2) ./ 2)).^2);

end