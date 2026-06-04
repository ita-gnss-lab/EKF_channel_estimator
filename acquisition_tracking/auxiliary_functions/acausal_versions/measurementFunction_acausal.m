function estimative = measurementFunction_acausal(x_k_k_1, configuration, q, localDelay)
%% Generate Autocorrelation Matrix
if nargin < 4 || isempty(localDelay)
    localDelay = 0;
end

delayError = x_k_k_1(1);
expectedTapCount = 2*q + 1;
if numel(x_k_k_1) == expectedTapCount + 1
    channelWeights = x_k_k_1(2:end);
    phaseScale = 1;
elseif numel(x_k_k_1) == expectedTapCount + 4
    channelWeights = x_k_k_1(5:end);
    phaseScale = exp(1j * x_k_k_1(2));
else
    error("measurementFunction_acausal:InvalidStateLength", ...
        "Expected state length %d or %d, got %d.", ...
        expectedTapCount + 1, expectedTapCount + 4, numel(x_k_k_1));
end
channelOrder = (numel(channelWeights) - 1)/2;

if isfield(configuration, "correlatorHalfSpan")
    correlatorHalfSpan = configuration.correlatorHalfSpan;
else
    correlatorHalfSpan = channelOrder;
end

shiftedCorrelationsMatrix = getShiftedCorrelations_acausal( ...
    delayError, correlatorHalfSpan, configuration, channelOrder, localDelay);

%% Get Estimated Measurement
estimative = phaseScale * ...
    sum(channelWeights.' .* shiftedCorrelationsMatrix, 2);

end
