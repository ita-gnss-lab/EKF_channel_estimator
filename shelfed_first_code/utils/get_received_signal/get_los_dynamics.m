function [los_phase, los_delay, time_vector] = get_los_dynamics( ...
        simulation_time, sampling_interval, doppler_profile, carrier_freq)
% get_los_dynamics
% Generates LOS carrier phase and corresponding delay time series using a
% Taylor expansion defined by `doppler_profile`.
%
% Syntax:
%   [los_phase, los_delay, t] = get_los_dynamics(Tsim, Ts, doppler_profile, fc)
%
% Inputs:
%   simulation_time   - Total duration [s], scalar > 0
%   sampling_interval - Sampling interval [s], scalar > 0
%   doppler_profile   - Row vector: [phase0, fd, fdr, coeff3, ...]
%                       phase0 in rad; fd in Hz; fdr in Hz/s; higher orders in Hz/s^(n-1)
%   carrier_freq      - Carrier frequency f_c [Hz], scalar > 0
%
% Outputs:
%   los_phase         - Column vector of carrier phase [rad]
%   los_delay         - Column vector of LOS delay [s]
%   time_vector       - Column vector of time instants [s]
%
% Notes:
%   Phase model:
%     phi(t) = phase0 + 2*pi * sum_{n=1}^{N-1} coeff_{n+1} * t^n / n!
%   Delay:
%     tau(t) = -phi(t) / (2*pi*fc)
%
% Author: Rodrigo de Lima Florindo
% ORCID: 0000-0003-0412-5583
% Email : rdlfresearch@gmail.com

% Input validation
validateattributes(simulation_time,  {'numeric'},{'scalar','real','positive','finite','nonnan'}, mfilename,'simulation_time');
validateattributes(sampling_interval, {'numeric'},{'scalar','real','positive','finite','nonnan'}, mfilename,'sampling_interval');
validateattributes(doppler_profile,   {'numeric'},{'row','real','finite','nonnan'},           mfilename,'doppler_profile');
validateattributes(carrier_freq,      {'numeric'},{'scalar','real','positive','finite','nonnan'}, mfilename,'carrier_freq');

if simulation_time < sampling_interval
    error('%s:simulationTimeSmallerThanSamplingInterval', mfilename, ...
          'simulation_time (%g) is smaller than sampling_interval (%g).', ...
           simulation_time, sampling_interval);
end

% Samples count (rounded if needed)
num_samples_exact   = simulation_time / sampling_interval;
num_samples_rounded = round(num_samples_exact);

if abs(num_samples_exact - num_samples_rounded) > eps
    warning('%s:NonIntegerRatio', mfilename, ...
        'simulation_time / sampling_interval not integer. Rounded from %.5g to %d samples.', ...
         num_samples_exact, num_samples_rounded);
end

% Time vector
time_vector = (1:num_samples_rounded).' * sampling_interval;

% Phase computation (vectorized Taylor series)
% Start with phase_0
los_phase = doppler_profile(1) * ones(size(time_vector));

% Add higher-order terms
% doppler_profile(2) is fd (1st order), so power starts at t^1 / 1!
orders = 2:length(doppler_profile);
if ~isempty(orders)
    for n = orders
        k = n - 1; % order in time exponent
        los_phase = los_phase + 2*pi * doppler_profile(n) .* (time_vector.^k) / factorial(k);
    end
end

% LOS delay from carrier phase
los_delay = -los_phase ./ (2*pi*carrier_freq);

end
