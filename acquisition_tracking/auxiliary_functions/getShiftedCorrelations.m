function shiftedCorrelations = getShiftedCorrelations(delayError, numberOfTaps, configuration)
correlationsDelay = delayError + ...
    1 / (configuration.samplingFrequency) * ...
    (-numberOfTaps : 1 : numberOfTaps);
shifts = 1 / (configuration.samplingFrequency) *...
    (0 : 1 : numberOfTaps);

shiftedDelayMatrix = correlationsDelay.' + shifts;

shiftedCorrelations = autocorrelation(shiftedDelayMatrix, configuration);

end