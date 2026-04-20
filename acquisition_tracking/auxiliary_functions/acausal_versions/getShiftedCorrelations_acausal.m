function shiftedCorrelations = getShiftedCorrelations_acausal(delayError, correlatorHalfSpan, configuration, channelOrder)
%GETSHIFTEDCORRELATIONS Build Φ_pp(ετ + (l - m)Ts) over correlator/tap indices.
%
%   delayError        : ετ term from the state vector.
%   correlatorHalfSpan: q parameter → correlator indices m ∈ [-q, q].
%   configuration     : receiver configuration structure.
%   channelOrder      : L parameter → channel taps l ∈ [0, L].
%
%   When channelOrder is omitted we assume L = q, matching the legacy
%   behaviour where the number of channel taps equalled the correlator
%   half-span.

    if nargin < 4 || isempty(channelOrder)
        channelOrder = correlatorHalfSpan;
    end

    Ts = 1 / configuration.samplingFrequency;

    m = (-correlatorHalfSpan:correlatorHalfSpan).';  % (2q+1)×1
    l = -channelOrder:channelOrder;                              % 1×(L+1)

    deltaMatrix = delayError + (l - m) * Ts;
    shiftedCorrelations = autocorrelation(deltaMatrix, configuration);
end
