
function y = rabifunc(x, c, d, B1, B0, Aa, Bb)
    if ~isempty(B0)
        y = (abs(c - d .* sin(atan2( B1./Aa, (x' - B0) )).^2.*sin(sqrt((B1./(Aa)).^2+(x' - B0).^2).*Bb).^2))';
    else
        y = [];
    end
end