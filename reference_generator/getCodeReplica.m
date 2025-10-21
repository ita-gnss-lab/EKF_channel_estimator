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
    maxDelay = ceil(abs(delaySamples)) + 2;
    delayHandler = dsp.VariableFractionalDelay( ...
        "InterpolationMethod","FIR", ...
        "MaximumDelay", maxDelay);

    signal = zeros(size(ref));
    reset(delayHandler);
    for idx = 1:length(ref)
        signal(idx) = delayHandler(ref(idx), delaySamples);
    end
end
