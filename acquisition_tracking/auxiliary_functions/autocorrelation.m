function [autocorrelation] = autocorrelation(delay, configuration)
%AUTOCORRELATION Summary of this function goes here
%   Detailed explanation goes here
epoch = configuration.totalChips / configuration.chippingFrequency;
signal = reference_signal(configuration, ...
                0)';
autocorrelation = zeros(size(delay, 1), size(delay, 2));
for i = 1:size(delay, 1)
    for j = 1:size(delay, 2)
    signalDelayed = reference_signal(configuration, ...
                    delay(i, j))';
    autocorrelation(i, j) = signal * signalDelayed.';
    end
end
end

