function [b1Out, b0Out, scores, errored] = radish(wnorm, mask, w, excludeB0)
%%  Parallel script for RADISH;
%   Inputs:
%       wnorm [ppm]: Normalised WASABI image (3D or 4D, x,y,(z),t)
%       mask [binary]: Mask image (2D or 3D, x,y,(z))
%       w: structure containing minimum peak prominence (mpp), scanner frequency in MHz (w0),
%           interpolation tolerance level (tol), and pulse duration in seconds (tp).
%   Outputs:
%       b1Out [uT]: Asbolute B1 map
%       b0Out [ppm]: B0 map
%       scores [au]: scores of the fit
%       errored [binary]: regions ignored due to error
%   Last updated March 2025 by Mara Quach (trinhq@student.unimelb.edu.au)


ROISize = size(wnorm);
if ndims(wnorm) == 4
    sz = ROISize(1) * ROISize(2) * ROISize(3);
    tt = 3;
elseif ndims(wnorm) == 3 
    sz = ROISize(1) * ROISize(2);
    tt = 2;
else
    disp("Please make sure input is either 3D (x,y,t) or 4D (x,y,z,t)")
end

% m1D = reshape(mask, sz, []);
m1D = mask(:);

if nargin < 4
    eb0 = 100 + zeros(sz,1);
else
    eb0 = excludeB0(:);
end

mpp = w.mpp;
w0 = w.w0;

wx = linspace(-w.max, w.max, w.noffsets);
wxi = linspace(-w.max, w.max, w.noffsets * w.ifactor);
sfactor = w.sfactor;
tp = w.tp; % make sure this is in s;
Aa = w0/w.gamma;
Bb = tp * w0 * pi;
Cc = (1/(tp * w0))^2;
lim = 2/tp/w0;
ps = wx(2) - wx(1);
if ps >= 0.0625 % equivalent to 49 offsets between -1.5 to 1.5
    ps = 0.0625;
end
mB1 = w.mB1; 
mB0 = w.mB0; % Should ensure that shimming is done 

% Making sure there are no NaNs
wnorm(isnan(wnorm)) = 0;
w1D = reshape(wnorm, sz,  []);

b1Out = zeros(ROISize(1:tt));
b0Out = zeros(ROISize(1:tt));
errored = zeros(ROISize(1:tt));
scores = zeros(ROISize(1:tt));

parfor jj = 1:sz % change to for if no parallel toolbox
    try
        if m1D(jj)
            vox = w1D(jj,:);
            sp = spaps(wx, vox, sfactor);
            [b1Out(jj), b0Out(jj), scores(jj)] = support_radish(vox', fnval(sp, wxi), wx, wxi, mB1, mB0, Aa, Bb, Cc, lim, mpp, ps, eb0(jj));
        end
    catch 
        errored(jj) = 1;
        continue
    end
end


end
