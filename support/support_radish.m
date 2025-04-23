function [ChosenB1, ChosenB0, score] = little_radish(vox, Y, wx, wxi, mB1, mB0, Aa, Bb, Cc, lim, mpp, ps, B0ignore)
    % Outputs:
    %           ChosenB1, ChosenB0: best B1, B0 candidates
    %           score: [au], double, score of the fit (lower = better)
    %           
    % Inputs:
    %           (name: [unit], datatype, description)
    %           vox: [ppm], 1xN double vector, where N = number of offsets 
    %           wx:  [ppm], 1xN double vector of offsets 
    %           wxi: [ppm], 1xN double vector, finely sampled offsets (suggested 10
    %                times original number of offsets) 
    %           w0:  [MHz], double, scanner frequency 
    %           mB1: [uT], 2x1 double vector, where mB1(1) is minimum accepted B1 and mB1(2)
    %                is maximum accepted B1 value
    %           mB0: [ppm], double, limit of B0 (e.g., mB0 = 0.8 -> -.8 <= B0 <= .8
    %                 accepted
    %           tp:  [s], double constant, saturation pulse duration
    %           Aa:  [T], double constant, field strength, or (w0/gamma).
    %           Bb:  [au], double constant, (tp * w0 * pi)
    %           Cc:  [au], double constant, ( (1/(tp * w0))^2 )
    %           lim: [ppm], double constant, peak distance limit, calculated by 2/tp/w0
    %           mpp: [au], double, minimum peak prominence. Recommended close
    %                 to 0
    %           ps:  [ppm], double, increment, difference between two consecutive data points
    %           sfact: [au], double, smoothing factor. Higher = smoother
    %           B0ignore: [ppm], double, B0 value to exclude
    
    
    % Find maxima/minima
    [~, px] = findpeaks(Y,'MinPeakProminence',mpp);
    [~, tx] = findpeaks(-Y,'MinPeakProminence',mpp);
        
    tx(wxi(tx) < -1.4 | wxi(tx) > 1.4) = [];
    
    nn = length(px);
    b0c1 = zeros(1, nn - 1); pd1 = zeros(1, nn - 1);
    b0c2 = zeros(1, nn - 2); pd2 = zeros(1, nn - 2);
    
    if nn > 4
        b0c3 = zeros(1, nn - 4); pd3 = zeros(1, nn - 4); lb1c3 = nn - 4;
    else
        b0c3 = zeros(1, 1); pd3 = zeros(1,1); lb1c3 = 1;
    end
    
    for ii = 1 : length(px) - 1
        cl = px(ii);
        cr = px(ii + 1); 
        pdi = wxi(cr) - wxi(cl);
        pd1(ii) = pdi;
    
        b0c1(ii) = wxi(cl) + pdi/2;
        
        if ii < length(px) - 1 
            cr = px(ii + 2);
            pdi = wxi(cr) - wxi(cl);
            pd2(ii) = pdi;
            b0c2(ii) = wxi(cl) + pdi/2;
            if ii < length(px) - 3 
                cr = px(ii + 4);
                pdi = wxi(cr) - wxi(cl);
                pd3(ii) = pdi;
                b0c3(ii) = wxi(cl) + pdi/2;
            end
        end 
    end
    
    b1c1 = findb1_l(pd1, Aa, Cc, lim);
    b1c2 = findb1_l(pd2, Aa, Cc, lim);
    b1c3 = findb1_l(pd3, Aa, Cc, lim);
    
    % clear cl cr pdi 
    
    lb1c1 = nn - 1;
    lb1c2 = nn - 2;
    
    ccurve = zeros(lb1c1 + lb1c2 + lb1c3 + 1, length(wx));
    
    ccurve(1:lb1c1, :) = rabifunc(wx, 1, 1, b1c1, b0c1, Aa, Bb);
    ccurve(lb1c1+1 : lb1c1+lb1c2, :) = rabifunc(wx, 1, 1.5, b1c2, b0c2, Aa, Bb);
    ccurve(lb1c1+lb1c2+1 : end - 1, :) = rabifunc(wx, 1, 2, b1c3, b0c3, Aa, Bb);
    
    b0c = [b0c1, b0c2, b0c3, 0];
    b1c = [b1c1, b1c2, b1c3, 0];
    pdc = [pd1, pd2, pd3, 0];
    
    % clear b0c1 b0c2 b0c3 b1c1 b1c2 b1c3 lb1c1 lb1c2 lb1c3 nn
    
    if length(tx) == 2 % cases where rB1 is very high (e.g., > 1.2)
        pxt = sort(reshape([tx px], [], 1));
        b0c(end) = (wxi(tx(2)) + wxi(tx(1))) / 2;
        pdc(end) = abs(b0c(end) - wxi(pxt(find(pxt == tx(2)) - 1)))*2;
        b1c(end) = findb1_l(pdc(end), Aa, Cc, lim);
        ccurve(end,:) = rabifunc(wx, 1, 1, b1c(end), b0c(end), Aa, Bb);
    end
    
    % clear tx pxt
    
    b0t = b1c~=0 & b0c > -mB0 & b0c > mB0(1) & b0c < mB0(2) & b0c~=B0ignore & b1c<=mB1(2) & b1c >=mB1(1);
    
    b0c = b0c(b0t);
    pdc = pdc(b0t);
    
    ccurve = ccurve(b0t,:);
    
    tcg = gradient(vox)';
    
    idx = pickle(tcg, ccurve);
    ChosenB0 = b0c(idx);
    pd = pdc(idx);
    
    if ismember(pd, pd1) 
        d = 1;
    elseif ismember(pd, pd2)
        d = 1.5;
    elseif ismember(pd, pd3)
        d = 2;
    else
        d = 1;
    end
    
    % clear pdc idx b0t pd1 pd2 pd3 ccurve
    
    if pd >= lim
        b1e = findb1_l(pd - ps*[1.5:.05:2.5], Aa, Cc, lim);
    else
        b1e = findb1_l(pd - ps*[.5:.05:1.5], Aa, Cc, lim);
    end
    
    b1e(b1e < mB1(1) | b1e > mB1(2)) = [];
    
    ccurvee = rabifunc(wx, 1, d, b1e, ChosenB0, Aa, Bb);
    
    [idx, score] = pickle(tcg, ccurvee);
    ChosenB1 = b1e(idx);
    
    end
