function [estimative, aux] = measurementFunction(x_k_k_1, configuration)
%% Generate Autocorrelation Matrix
delayError = x_k_k_1(1);
channelWeights = x_k_k_1(5:end);
channelOrder = numel(channelWeights) - 1;
dopplerError = x_k_k_1(3);
dopplerRateError = x_k_k_1(4);

if isfield(configuration, "correlatorHalfSpan")
    correlatorHalfSpan = configuration.correlatorHalfSpan;
else
    correlatorHalfSpan = channelOrder;
end

shiftedCorrelationsMatrix = getShiftedCorrelations( ...
    delayError, correlatorHalfSpan, configuration, channelOrder);

Ts = 1 / configuration.samplingFrequency;
tapIndices = (0:channelOrder).';
tapTimes = delayError + tapIndices * Ts;

dopplerPhase = exp(1j * (dopplerError * tapTimes + ...
    0.5 * dopplerRateError * tapTimes.^2));
effectiveWeights = channelWeights .* dopplerPhase;

%% Get Estimated Measurement
totalPhaseError = x_k_k_1(2);
estimative = exp(1j * totalPhaseError) * ...
    sum(effectiveWeights.' .* shiftedCorrelationsMatrix, 2); 

if nargout > 1
    aux.shiftedCorrelations = shiftedCorrelationsMatrix;
    aux.tapTimes = tapTimes;
    aux.dopplerPhase = dopplerPhase;
    aux.effectiveWeights = effectiveWeights;
    aux.totalPhaseError = totalPhaseError;
end

end
