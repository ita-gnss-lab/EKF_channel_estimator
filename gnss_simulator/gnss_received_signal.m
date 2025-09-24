function [received_signal, time] = gnss_received_signal(configuration, numberOfEpochs)
%   GNSS_RECEIVED_SIGNAL
%   Detailed explanation goes here

N = numberOfEpochs * configuration.totalChips * ...
            (configuration.samplingFrequency / configuration.chippingFrequency);
% Time vector
time = (0 :  N - 1) / configuration.samplingFrequency;

%% TRANSMITTED PILOT SIGNAL
% Gets the code in bits; 0 and 1
PRN_code = gnssCACode(configuration.satellite, "GPS");
% Remaps to +1 and -1
ranging_code = double(2*PRN_code - 1);
% Samples it
current_chip = floor(mod((time.* configuration.chippingFrequency), numel(PRN_code)) + 1);
sampled_code = ranging_code(current_chip);  

%% LOS PHASE AND DELAY
% Computes LOS transmission phase and delay
[LOS_phase, LOS_delay] = get_LOS_dynamics(time, configuration.dopplerProfile, configuration.carrierFrequency);
% Applies the delay
% delay_in_samples = LOS_delay * configuration.samplingFrequency;
% delay_handler = dsp.VariableFractionalDelay("InterpolationMethod","Linear", 'MaximumDelay',9999);
% delayed_code = delay_handler(sampled_code, delay_in_samples);
delay_in_samples = floor(LOS_delay * configuration.samplingFrequency);
delayed_code = circshift(sampled_code, delay_in_samples);
% Applies the phase
received_signal = delayed_code .* exp(1j*LOS_phase);

%% THERMAL NOISE
% noise = get_simple_thermal_noise(length(time), 1 / configuration.samplingFrequency, configuration.carrierToNoiseDensityRatio);
% received_signal = distorted_code + noise; 

end

