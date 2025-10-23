function signal = getCodeReplica(configuration, delaySeconds)
    Fs = configuration.samplingFrequency;
    Fc = configuration.chippingFrequency;
    K  = configuration.totalChips;

    code = gnssCACode(configuration.satellite, "GPS");
    codeBPSK = 2*double(code) - 1;

    Nfloat = K * (Fs / Fc);
    N = round(Nfloat);
    if abs(N - Nfloat) > 1e-9
        warning('Non-integer samples/epoch (%.9f). Using N=%d.', Nfloat, N);
    end

    n = 0:(N-1);
    idx0 = floor(mod(n * (Fc/Fs), K)) + 1;
    ref  = codeBPSK(idx0);

    delaySamples = delaySeconds * Fs;

    % Wrap the evaluation grid and use interp1 for the fractional delay.
    % Append the first sample to preserve periodicity during interpolation.
    queryPoints = mod(((0:N-1).' - delaySamples), N);
    refExtended = [ref(:); ref(1)];

    signal = interp1(0:N, refExtended, queryPoints, "linear");
    signal = signal(:).';
end
