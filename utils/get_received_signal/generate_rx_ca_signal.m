function [rx_sig, ca_bb, ca_bb_del, los_phase, los_delay, time] = ...
    generate_rx_ca_signal(prn, C_over_N0_dBHz, doppler_profile, fc, simulation_time, Ts)
% generate_rx_ca_signal
% 
% Build a noisy complex baseband GPS L1 C/A signal (single PRN).
%
% Inputs
%   prn               : PRN number (1..32)
%   C_over_N0_dBHz    : C/N0 in dB-Hz
%   doppler_profile   : [phi0, fd, fdr, ...] for get_los_dynamics
%   fc                : carrier frequency (e.g. 1.57542e9 Hz)
%   simulation_time   : total duration [s]
%   Ts                : sampling interval [s]  (Fs = 1/Ts must be integer multiple of chip rate/1.023e6)
%
% Outputs
%   rx_sig            : complex baseband received signal (C/A * e^{j phi_LOS} + noise)
%   ca_bb             : upsampled C/A code (real, +/-1)
%   los_phase         : LOS carrier phase (rad)
%   los_delay         : LOS delay (s)
%   time              : time vector [s]
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

%% Parameters
Rc = 1.023e6;                 % C/A chip rate [Hz]
Fs = 1/Ts;                    % sample rate [Hz]
samples_per_chip = Fs/Rc;

if abs(samples_per_chip - round(samples_per_chip)) > 1e-12
    error('Fs must be an integer multiple of 1.023e6. Got %.6f samples/chip.', samples_per_chip);
end
samples_per_chip = round(samples_per_chip);

%% phase / delay / time
[los_phase, los_delay, time] = get_los_dynamics(simulation_time, Ts, doppler_profile, fc);

%% Generate upsampled C/A code (pilot, no data)
chips = gnssCACode(prn, "GPS"); % 0/1
chips = double(2*chips - 1);    % +/-1
% Repeat to cover whole simulation
N     = numel(time);
idx   = mod(floor((0:N-1).'/samples_per_chip), numel(chips)) + 1;
ca_bb = chips(idx);                  % column

% Desired delay in samples (vector)
d = los_delay * Fs;                 % total delay [samples], N×1

vfd = dsp.VariableFractionalDelay("InterpolationMethod","Linear", 'MaximumDelay',9999);

ca_bb_del = vfd(ca_bb, d);            % same length as ca_bb

%% LOS modulation
sig_los = ca_bb_del .* exp(1j*los_phase);  % complex baseband

%% Thermal noise
noise = get_thermal_noise(simulation_time, Ts, C_over_N0_dBHz);

%% Received signal
rx_sig = sig_los + noise;

end
