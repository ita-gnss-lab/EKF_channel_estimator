function mape = meanAbsolutePercentageError(yTrue, yPred)
%MEANABSOLUTEPERCENTAGEERROR Uniform mean absolute percentage error.

denominator = max(abs(yTrue), eps);
relativeError = abs(yPred - yTrue) ./ denominator;
mape = 100 * mean(relativeError(:));

end
