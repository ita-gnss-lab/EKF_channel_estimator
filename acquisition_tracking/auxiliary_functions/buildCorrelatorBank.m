function bank = buildCorrelatorBank(cfg, delaySeconds, q)
% Integer-sample correlator bank: each row m is the prompt code
% shifted by m samples, m ∈ [-q..q].  (Δτ = Ts)

    Fs = cfg.samplingFrequency;
    Fc = cfg.chippingFrequency;
    K  = cfg.totalChips;                 % 1023 for GPS C/A

    % --- PRN in ±1
    code01 = gnssCACode(cfg.satellite, "GPS");
    codePM = 2*double(code01) - 1;       % length K

    % --- samples per epoch (round once)
    Nfloat = K * (Fs / Fc);
    N = round(Nfloat);
    if abs(N - Nfloat) > 1e-9
        warning('Non-integer samples/epoch (%.9f). Using N=%d.', Nfloat, N);
    end

    % --- sample the zero-phase code to N samples
    n = 0:(N-1);
    idx0 = floor(mod(n * (Fc/Fs), K)) + 1;  % 1..K
    ref  = codePM(idx0);                    % 1×N zero-phase sampled code

    % --- make prompt by quantized (integer-sample) a priori delay
    k0 = round(delaySeconds * Fs);          % samples
    prompt = circshift(ref, k0);            % 1×N

    % Build bank by integer shifts m ∈ [-q .. q]
    bank = zeros(2*q + 1, N);   % preallocate (rows = shifts, columns = samples)
    
    row = 1;
    for m = -q:q
        % Shift the prompt by m samples (positive m = shift right)
        bank(row, :) = circshift(prompt, m);
        row = row + 1;
    end
end
