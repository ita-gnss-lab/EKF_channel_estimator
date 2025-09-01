function [shiftedDerivatives] = getShiftedDerivatives(delayError, numberOfTaps, configuration)
%GETSHIFTEDDERIVATIVES Summary of this function goes here
%   Detailed explanation goes here
samplingPeriod = 1 / configuration.samplingFrequency;
correlationsDelay = delayError + samplingPeriod * (-numberOfTaps : 1 : numberOfTaps);
shifts = samplingPeriod * (0 : 1 : numberOfTaps);

shiftedDelayMatrix = correlationsDelay' + shifts;
shiftedDerivatives = zeros(size(shiftedDelayMatrix, 1), size(shiftedDelayMatrix, 2));
for i = 1:size(shiftedDelayMatrix, 1)
    for j = 1:size(shiftedDelayMatrix, 2)
    shiftedDerivatives(i, j) = autocorrelationDerivative(shiftedDelayMatrix(i, j), ...
        configuration.chippingFrequency, configuration.samplingFrequency);
    end
end
end

