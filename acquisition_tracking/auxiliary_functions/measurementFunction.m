function estimative = measurementFunction(x_k_k_1, configuration)
%% Generate Autocorrelation Matrix
numberOfTaps = length(x_k_k_1(6:end));
delayError = x_k_k_1(1);
shiftedCorrelationsMatrix = ...
    getShiftedCorrelations(delayError, numberOfTaps, configuration);

%% Get Estimated Measurement
totalPhaseError = x_k_k_1(2);
channelWeights = x_k_k_1(5:end);
estimative = exp(1j * totalPhaseError) * ...
    sum(channelWeights.' .* shiftedCorrelationsMatrix, 2); 

end