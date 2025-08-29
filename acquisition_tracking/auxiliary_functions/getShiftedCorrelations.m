function shiftedCorrelations = getShiftedCorrelations(delayError, numberOfTaps, samplingPeriod)

correlationsDelay = delayError + samplingPeriod * (-numberOfTaps : 1 : numberOfTaps);
shifts = samplingPeriod * (0 : 1 : numberOfTaps);

shiftedDelayMatrix = correlationsDelay' + shifts;

shiftedCorrelations = autocorrelation(shiftedDelayMatrix);

end