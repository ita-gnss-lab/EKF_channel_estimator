function bank = buildCorrelatorBank(configuration, delaySeconds, q)
% Integer-sample correlator bank: each row m is the prompt code
% shifted by m samples, m ∈ [-q..q].  (Δτ = Ts)

    signal = getCodeReplica(configuration, delaySeconds);

    row = 1;
    for m = -q:q
        % Shift the prompt by m samples (positive m = shift right)
        bank(row, :) = circshift(signal, m);
        row = row + 1;
    end
end
