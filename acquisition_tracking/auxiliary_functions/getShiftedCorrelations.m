function shiftedCorrelations = getShiftedCorrelations(delayError, numberOfTaps, configuration)
correlationsDelay = delayError + ...
    1 / (configuration.chippingFrequency * numberOfTaps) * ...
    (-numberOfTaps : 1 : numberOfTaps);
shifts = 1 / (configuration.chippingFrequency * numberOfTaps) *...
    (0 : 1 : numberOfTaps);

shiftedDelayMatrix = correlationsDelay' + shifts;

shiftedCorrelations = autocorrelation(shiftedDelayMatrix, configuration);

end