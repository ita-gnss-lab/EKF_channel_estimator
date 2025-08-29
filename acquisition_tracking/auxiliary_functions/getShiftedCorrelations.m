function shiftedCorrelations = getShiftedCorrelations(delayError, numberOfTaps, configuration)
samplingPeriod = 1 / configuration.samplingFrequency;
correlationsDelay = delayError + samplingPeriod * (-numberOfTaps : 1 : numberOfTaps);
shifts = samplingPeriod * (0 : 1 : numberOfTaps);

shiftedDelayMatrix = correlationsDelay' + shifts;

shiftedCorrelations = autocorrelation(shiftedDelayMatrix, configuration);

end