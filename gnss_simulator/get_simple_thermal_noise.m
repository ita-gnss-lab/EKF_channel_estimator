function thermal_noise = get_simple_thermal_noise(simulation_time, Ts, C_over_N0_dBHz)
% get_thermal_noise
% Generates additive white Gaussian noise (AWGN) to simulate thermal noise 
% in a receiver based on specified parameters.
%
% Syntax:
%   thermal_noise = get_thermal_noise(simulation_time, Ts, C_over_N0_dBHz)
%
% Description:
%   This function generates additive white Gaussian noise (AWGN) to model 
%   thermal noise in a receiver. The noise is simulated as a wide-sense 
%   stationary complex Gaussian random process. Its variance is determined 
%   based on the carrier-to-noise density ratio (C/N₀), and correlation 
%   sampling interval.
%
% Inputs:
%   simulation_time - Total duration of the simulation (seconds).
%   Ts - Sampling interval (seconds).
%   C_over_N0_dBHz  - Carrier-to-noise density ratio (C/N₀) in dB-Hz.
%
% Outputs:
%   thermal_noise   - Complex Gaussian noise time series with a variance 
%                     derived from the specified parameters.
%
% Notes:
%   - Consider the complex baseband signal sampled at Fs = 1/Ts = 2B (Hz). 
%   If the received signal power is normalized to 1, the continuous-time 
%   (per-Hz) noise PSD is 1 / (c_over_n0). Spreading this over the sampled 
%   bandwidth gives the per-sample complex noise variance:
%
%       σ_η_d² = 1 / (c_over_n0 * Ts) (1)
%
% Example:
%   % Generate thermal noise for:
%   % - simulation_time = 100 seconds
%   % - Ts = 1/4e6 seconds
%   % - C_N0 = 40 dB-Hz
%   thermal_noise = get_thermal_noise(simulation_time, Ts, C_N0);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

% % Input validation
% validateattributes(simulation_time, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'},'get_thermal_noise', 'simulation_time');
% validateattributes(Ts, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_thermal_noise', 'Ts');
% validateattributes(C_over_N0_dBHz, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_thermal_noise', 'C_over_N0_dBHz');
% 
% if simulation_time < Ts
%     error('get_thermal_noise:simulationTimeSmallerThanSamplingInterval', ...
%         'The inputed value of `simulation_time` was %g, which is smaller than the value of the `Ts`, %g', simulation_time, Ts)
% end
% 
% % Check if simulation_time / Ts is an integer
% num_samples_exact = simulation_time / Ts;
% num_samples_rounded = round(num_samples_exact);
% 
% if abs(num_samples_exact - num_samples_rounded) > eps
%     % Issue a warning if rounding was needed
%     warning('get_thermal_noise:NonIntegerRatioSamples', ...
%             'simulation_time / Ts is not an integer. The number of samples was rounded from %.5g to %d.', ...
%             num_samples_exact, num_samples_rounded);
% end
num_samples_rounded = simulation_time; 
% Convert CN0 from dB-Hz to linear scale
c_over_n0_linear = 10^(C_over_N0_dBHz / 10);

% Compute the noise variance
sigma2_eta_d = 1 / (Ts * c_over_n0_linear);

% Generate complex Gaussian noise
real_part = randn(num_samples_rounded, 1) * sqrt(sigma2_eta_d / 2);
imag_part = randn(num_samples_rounded, 1) * sqrt(sigma2_eta_d / 2);

% Combine real and imaginary parts
thermal_noise = real_part + 1j * imag_part;

end
