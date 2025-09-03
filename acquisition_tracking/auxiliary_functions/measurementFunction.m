function estimative = measurementFunction(stateAPriori, configuration)
%% Generate Autocorrelation Matrix
numberOfTaps = length(stateAPriori(6:end));
delayError = stateAPriori(1);
shiftedCorrelationsMatrix = ...
    getShiftedCorrelations(delayError, numberOfTaps, configuration);

%% Get Estimated Measurement
totalPhaseError = stateAPriori(2);
channelWeights = stateAPriori(5:end);
estimative = exp(1j * totalPhaseError) * ...
    sum(channelWeights.' .* shiftedCorrelationsMatrix, 2); 

end