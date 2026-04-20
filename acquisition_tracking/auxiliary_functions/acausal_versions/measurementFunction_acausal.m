function estimative = measurementFunction_acausal(x_k_k_1, configuration, q)
%% Generate Autocorrelation Matrix
delayError = x_k_k_1(1);
channelWeights = x_k_k_1(5:end);
channelOrder = (numel(channelWeights) - 1)/2;
tapOrder = -q:q;
tapPhase = (-(2*pi*configuration.carrierFrequency)/configuration.samplingFrequency)*tapOrder;
channelWeights = channelWeights .* exp(1j * tapPhase');

if isfield(configuration, "correlatorHalfSpan")
    correlatorHalfSpan = configuration.correlatorHalfSpan;
else
    correlatorHalfSpan = channelOrder;
end

shiftedCorrelationsMatrix = getShiftedCorrelations_acausal( ...
    delayError, correlatorHalfSpan, configuration, channelOrder);

% epoch = configuration.totalChips / configuration.chippingFrequency;
% samplesTotal = epoch*configuration.samplingFrequency;
% n = 0:(samplesTotal-1);
% time = n/configuration.samplingFrequency;

% dopplerProfile = [x_k_k_1(2), ...
%     x_k_k_1(3)/(2*pi),...
%     x_k_k_1(4)/(2*pi)];
% totalPhaseError = get_LOS_dynamics(time, dopplerProfile, configuration.carrierFrequency);

%% Get Estimated Measurement
totalPhaseError = x_k_k_1(2);
estimative = exp(1j * totalPhaseError) * ...
    sum(channelWeights.' .* shiftedCorrelationsMatrix, 2); 

end
