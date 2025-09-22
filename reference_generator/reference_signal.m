function signal=reference_signal(configuration, delay)
% Time vector
numberOfEpochs = 1;
samples = numberOfEpochs * configuration.totalChips * ...
            (configuration.samplingFrequency / configuration.chippingFrequency);

% Time vector
time = 1 / configuration.samplingFrequency * ...
        (0 :  samples);

%% TRANSMITTED PILOT SIGNAL
% Gets the code in bits; 0 and 1
PRN_code = gnssCACode(configuration.satellite, "GPS");
% Remaps to +1 and -1
ranging_code = double(2*PRN_code - 1);
% Samples it
current_chip = mod(floor(time.* configuration.chippingFrequency), numel(PRN_code)) + 1;
sampled_code = ranging_code(current_chip);  
% Applies the delay
delay_in_samples = floor(delay * configuration.samplingFrequency);
signal = circshift(sampled_code, delay_in_samples);

