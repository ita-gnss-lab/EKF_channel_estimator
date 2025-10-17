function signal = getCodeReplica(configuration, delaySeconds)
    Fs = configuration.samplingFrequency;
    Fc = configuration.chippingFrequency;
    K  = configuration.totalChips;                 % 1023 for GPS C/A
    
    % --- PRN in ±1
    code = gnssCACode(configuration.satellite, "GPS");
    codeBPSK = 2*double(code) - 1;       % length K
    
    % --- samples per epoch (round once)
    Nfloat = K * (Fs / Fc);
    N = round(Nfloat);
    if abs(N - Nfloat) > 1e-9
        warning('Non-integer samples/epoch (%.9f). Using N=%d.', Nfloat, N);
    end
    
    % --- sample the zero-phase code to N samples
    n = 0:(N-1);
    idx0 = floor(mod(n * (Fc/Fs), K)) + 1;  % 1..K
    ref  = codeBPSK(idx0);                    % 1×N zero-phase sampled code
    
    % --- make prompt by quantized (integer-sample) a priori delay
    k0 = round(delaySeconds * Fs);          % samples
    signal = circshift(ref, k0);            % 1×N
end
