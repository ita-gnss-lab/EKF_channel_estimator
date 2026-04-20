function bank = buildCorrelatorBank(configuration, delaySeconds, q)
% Integer-sample correlator bank: each row m is the prompt code
% shifted by m samples, m ∈ [-q..q].  (Δτ = Ts)

    signal = getCodeReplica(configuration, delaySeconds);
    Fs = configuration.samplingFrequency;
    Fc = configuration.chippingFrequency;
    K  = configuration.totalChips;
    Nfloat = K * (Fs / Fc);
    N = round(Nfloat);
    row = 1;
    bank = zeros(2*q+1, N);
    for m = -q:q
        % Shift the prompt by m samples (positive m = shift right)
        bank(row, :) = circshift(signal, m);
        row = row + 1;
    end
end
