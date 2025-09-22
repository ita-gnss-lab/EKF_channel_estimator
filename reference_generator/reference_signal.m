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



% % Computes LOS transmission phase and delay
% doppler = configuration.dopplerProfile(2:end);
% [LOS_phase, LOS_delay] = get_LOS_dynamics(time, [delay doppler], configuration.carrierFrequency);
% Applies the delay
delay_in_samples = delay * configuration.samplingFrequency;
delay_handler = dsp.VariableFractionalDelay("InterpolationMethod","Linear", 'MaximumDelay',9999);
signal_delayed = delay_handler(sampled_code, delay_in_samples);
signal_circ_shifted = circshift(signal_delayed, round(delay_in_samples));
signal_delayed(1:round(delay_in_samples)) = signal_circ_shifted(1:round(delay_in_samples));
signal = signal_delayed;


% %Applies the delay
% delay_in_samples = floor(delay * configuration.samplingFrequency);
% signal = circshift(sampled_code, delay_in_samples);

