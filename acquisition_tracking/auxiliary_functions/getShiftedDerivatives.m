function [shiftedDerivatives] = getShiftedDerivatives(delayError, numberOfTaps, samplingPeriod)
%GETSHIFTEDDERIVATIVES Summary of this function goes here
%   Detailed explanation goes here
correlationsDelay = delayError + samplingPeriod * (-numberOfTaps : 1 : numberOfTaps);
shifts = samplingPeriod * (0 : 1 : numberOfTaps);

shiftedDelayMatrix = correlationsDelay' + shifts;

shiftedDerivatives = autocorrelationDerivative(shiftedDelayMatrix, ...
    configuration.chippingFrequency, configuration.samplingFrequency);
end

