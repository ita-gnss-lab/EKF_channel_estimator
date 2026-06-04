function shiftedCorrelations = getShiftedCorrelationsFromBank( ...
    correlatorBank, delay, configuration, channelTapDelays)
%GETSHIFTEDCORRELATIONSFROMBANK Correlate bank rows with causal tap replicas.

samplesTotal = size(correlatorBank, 2);
channelReplicas = zeros(numel(channelTapDelays), samplesTotal);
for channelIndex = 1:numel(channelTapDelays)
    channelReplicas(channelIndex, :) = getCodeReplica( ...
        configuration, delay + channelTapDelays(channelIndex)).';
end

shiftedCorrelations = correlatorBank * channelReplicas.';

end
