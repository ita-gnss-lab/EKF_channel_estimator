function bank = buildCorrelatorBank(configuration, delaySeconds, q)
% Integer-sample correlator bank: row m corresponds to the prompt code
% shifted by -m samples, with m ∈ {q, ..., -q}. (Delta tau multiples of Ts)

    signal = getCodeReplica(configuration, delaySeconds);

    mValues = q:-1:-q;
    row = 1;
    for m = mValues
        % Shift replica by -m so Delta_{l,m} = eps_tau + (l + m)Ts in the model
        % NOTE: For circshift, the second argument is the delay in samples.
        % Therefore, +1 in the input generates a -1 sample delay. This is
        % way we need to use -m.
        bank(row, :) = circshift(signal, -m);
        row = row + 1;
    end
end
