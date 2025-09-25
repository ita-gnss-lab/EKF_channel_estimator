function [receivedSignal, time] = gnss_received_signal(configuration, numberOfEpochs)
%   GNSS_RECEIVED_SIGNAL
%   Detailed explanation goes here

N = numberOfEpochs * configuration.totalChips * ...
            (configuration.samplingFrequency / configuration.chippingFrequency);
% Time vector
time = (0 :  N - 1) / configuration.samplingFrequency;

%% TRANSMITTED PILOT SIGNAL
% Gets the code in bits; 0 and 1
PRNCode = gnssCACode(configuration.satellite, "GPS");
% Remaps to +1 and -1
rangingCode = double(2*PRNCode - 1);
% Samples it
currentChip = floor(mod((time * configuration.chippingFrequency), numel(PRNCode)) + 1);
sampledCode = rangingCode(currentChip);  

%% LOS PHASE AND DELAY
% Computes LOS transmission phase and delay
[LOSPhase, LOSDelay] = get_LOS_dynamics(time, configuration.dopplerProfile, configuration.carrierFrequency);
% Applies the delay
% delay_in_samples = LOS_delay * configuration.samplingFrequency;
% delay_handler = dsp.VariableFractionalDelay("InterpolationMethod","Linear", 'MaximumDelay',9999);
% delayed_code = delay_handler(sampled_code, delay_in_samples);
delayInSamples = round(LOSDelay * configuration.samplingFrequency);
delayedCode = circshift(sampledCode, delayInSamples);
% Applies the phase
receivedSignal = delayedCode .* exp(1j*LOSPhase);

%% THERMAL NOISE
% noise = get_simple_thermal_noise(length(time), 1 / configuration.samplingFrequency, configuration.carrierToNoiseDensityRatio);
% received_signal = distorted_code + noise; 

end

