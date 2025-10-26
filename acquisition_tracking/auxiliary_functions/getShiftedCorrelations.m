function shiftedCorrelations = getShiftedCorrelations(delayError, q, configuration)
%GETSHIFTEDCORRELATIONS Build Φ_pp(ετ + (l - m)Ts) over correlator/tap indices.
%
%   delayError: ετ term from the state vector.
%   correlatorHalfSpan: q parameter - correlator indices m ∈ [-q, q].
%   configuration: receiver configuration structure.

    Ts = 1 / configuration.samplingFrequency;

    m = (-q:q).';  % (2q+1)×1
    l = 0:q;                              % 1×(L+1)

    deltaMatrix = delayError + (l - m) * Ts;
    shiftedCorrelations = autocorrelation(deltaMatrix, configuration);
end
