function [jacobian] = delayJacobianFunction(stateAPriori, samplingPeriod)
%DELAYJACOBIANFUNCTION Summary of this function goes here
%   Detailed explanation goes here
%% Generate Autocorrelation Matrix
numberOfTaps = length(stateAPriori(5:end));
delayError = stateAPriori(1);
shiftedDerivativesMatrix = ...
    getShiftedDerivatives(delayError, numberOfTaps, samplingPeriod);

%% Get Estimated Measurement
totalPhaseError = stateAPriori(2);
channelWeights = stateAPriori(5:end);
jacobian = exp(totalPhaseError) * ...
    sum(channelWeights .* shiftedDerivativesMatrix); 
end
