function w = setup_radish
    w.sfactor = 3e-4; % smoothing factor
    w.mpp = 0; % minimum peak prominence [au]
    w.gamma = 42.576; % Gyromagnetic ratio [MHz/T]
    w.ifactor = 1000; % interpolation factor
    w.mB1 = [0 5.5]; % B1 range [uT]
    w.mB0 = [-.6 .6]; % B0 range [ppm]

    w.tp = nan; % s
    w.w0 = nan; % MHz
    w.max = nan; % maximum offset value [ppm]
    w.noffsets = nan; % number of offsets
    display("Struct initialised. Please set these values accordingly:" + newline + " > tp = pulse duration [s]" + newline + " > w0: Scanner frequency [MHz]" + newline + " > max: maximum offset value [ppm]" + newline + " > noffsets: number of offsets");
end