function [autocorrelation] = autocorrelation(delay, configuration)
%AUTOCORRELATION Summary of this function goes here
%   Detailed explanation goes here
signal = getCodeReplica(configuration, 0)';
autocorrelation = zeros(size(delay, 1), size(delay, 2));
for i = 1:size(delay, 1)
    for j = 1:size(delay, 2)
    signalDelayed = getCodeReplica(configuration, delay(i, j))';
    autocorrelation(i, j) = signal * signalDelayed.';
    end
end
end

