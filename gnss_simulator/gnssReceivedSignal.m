function [receivedSignal, time, LOSPhase, LOSDelay] = gnssReceivedSignal(configuration, numberOfEpochs)
%   GNSS_RECEIVED_SIGNAL
%   Detailed explanation goes here

N = numberOfEpochs * configuration.totalChips * ...
            (configuration.samplingFrequency / configuration.chippingFrequency);
% Time vector
n = 0:(N-1);
time = n / configuration.samplingFrequency;

%% TRANSMITTED PILOT SIGNAL
% Gets the code in bits; 0 and 1
PRNCode = gnssCACode(configuration.satellite, "GPS");
% Remaps to +1 and -1
rangingCode = double(2*PRNCode - 1);
% Samples it
% NOTE(Rodrigo): It was not okay to use time *% configuration.chippingFrequency. 
% That is because time was a float, and it was causing some rounding errors
% sometimes. This was the cause of the incorrect triangle shape that we 
% were observing in the plots.
currentChip = floor(mod(n * (configuration.chippingFrequency / configuration.samplingFrequency), numel(PRNCode)) + 1);
sampledCode = rangingCode(currentChip);  

%% LOS PHASE AND DELAY
% Computes LOS transmission phase and delay
[LOSPhase, LOSDelay] = get_LOS_dynamics(time, configuration.dopplerProfile, configuration.carrierFrequency);

delayInSamples = LOSDelay * configuration.samplingFrequency;
maxDelay = ceil(max(delayInSamples)) + 2;
delayHandler = dsp.VariableFractionalDelay( ...
    "InterpolationMethod","FIR", ...
    "MaximumDelay", maxDelay);

delayedCode = zeros(size(sampledCode));
reset(delayHandler);
for idx = 1:length(sampledCode)
    delayedCode(idx) = delayHandler(sampledCode(idx), delayInSamples(idx));
end
% Applies the phase
distorted_code = delayedCode .* exp(1j*LOSPhase);

%% THERMAL NOISE
if configuration.addNoise 
    noise = get_simple_thermal_noise(length(time), 1 / configuration.samplingFrequency, configuration.carrierToNoiseDensityRatio);
    receivedSignal = distorted_code + noise; 
else
    receivedSignal = distorted_code;
end

end
