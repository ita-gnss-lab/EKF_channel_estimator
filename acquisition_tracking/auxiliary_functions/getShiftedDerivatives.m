function [shiftedDerivatives] = getShiftedDerivatives(delayError, numberOfTaps, configuration)
%GETSHIFTEDDERIVATIVES Summary of this function goes here
%   Detailed explanation goes here
correlationsDelay = delayError + ...
    1 / (configuration.chippingFrequency * numberOfTaps) * ...
    (-numberOfTaps : 1 : numberOfTaps);
shifts = 1 / (configuration.chippingFrequency * numberOfTaps) *...
    (0 : 1 : numberOfTaps);

shiftedDelayMatrix = correlationsDelay.' + shifts;
shiftedDerivatives = zeros(size(shiftedDelayMatrix, 1), size(shiftedDelayMatrix, 2));
for i = 1:size(shiftedDelayMatrix, 1)
    for j = 1:size(shiftedDelayMatrix, 2)
    shiftedDerivatives(i, j) = autocorrelationDerivative(shiftedDelayMatrix(i, j), ...
        configuration.chippingFrequency, configuration.samplingFrequency);
    end
end
end

