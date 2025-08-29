function estimative = measurementFunction(stateAPriori, samplingPeriod)
%% Generate Autocorrelation Matrix
numberOfTaps = length(stateAPriori(5:end));
delayError = stateAPriori(1);
shiftedCorrelationsMatrix = ...
    getShiftedCorrelations(delayError, numberOfTaps, samplingPeriod);

%% Get Estimated Measurement
totalPhaseError = stateAPriori(2);
channelWeights = stateAPriori(5:end);
estimative = exp(totalPhaseError) * ...
    sum(channelWeights .* shiftedCorrelationsMatrix); 

end