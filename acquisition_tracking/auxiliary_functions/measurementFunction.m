function estimative = measurementFunction(x_k_k_1, configuration)
%% Generate Autocorrelation Matrix
delayError = x_k_k_1(1);
channelWeights = x_k_k_1(5:end);
channelOrder = numel(channelWeights) - 1;

if isfield(configuration, "correlatorHalfSpan")
    correlatorHalfSpan = configuration.correlatorHalfSpan;
else
    correlatorHalfSpan = channelOrder;
end

shiftedCorrelationsMatrix = getShiftedCorrelations( ...
    delayError, correlatorHalfSpan, configuration, channelOrder);

%% Get Estimated Measurement
totalPhaseError = x_k_k_1(2);
estimative = exp(1j * totalPhaseError) * ...
    sum(channelWeights.' .* shiftedCorrelationsMatrix, 2); 

end
