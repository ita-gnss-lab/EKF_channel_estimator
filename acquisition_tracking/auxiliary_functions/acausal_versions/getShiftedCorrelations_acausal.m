function shiftedCorrelations = getShiftedCorrelations_acausal(delayError, correlatorHalfSpan, configuration, channelOrder, localDelay)
%GETSHIFTEDCORRELATIONS Build Φ_pp(ετ + (l - m)Ts) over correlator/tap indices.
%
%   delayError        : ετ term from the state vector.
%   correlatorHalfSpan: q parameter → correlator indices m ∈ [-q, q].
%   configuration     : receiver configuration structure.
%   channelOrder      : L parameter → channel taps l ∈ [-L, L].
%   localDelay        : absolute delay of the local prompt replica.
%
%   When channelOrder is omitted we assume L = q, matching the legacy
%   behaviour where the number of channel taps equalled the correlator
%   half-span.

    if nargin < 4 || isempty(channelOrder)
        channelOrder = correlatorHalfSpan;
    end
    if nargin < 5 || isempty(localDelay)
        localDelay = 0;
    end

    Ts = 1 / configuration.samplingFrequency;
    samplesPerEpoch = round(configuration.totalChips * ...
        (configuration.samplingFrequency / configuration.chippingFrequency));

    m = (-correlatorHalfSpan:correlatorHalfSpan).';
    l = -channelOrder:channelOrder;

    correlatorReplicas = zeros(numel(m), samplesPerEpoch);
    for row = 1:numel(m)
        correlatorReplicas(row, :) = ...
            getCodeReplica(configuration, localDelay + m(row)*Ts).';
    end

    channelReplicas = zeros(numel(l), samplesPerEpoch);
    for col = 1:numel(l)
        channelReplicas(col, :) = ...
            getCodeReplica(configuration, localDelay + delayError + l(col)*Ts).';
    end

    shiftedCorrelations = correlatorReplicas * channelReplicas.';
end
