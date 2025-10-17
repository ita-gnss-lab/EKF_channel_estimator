configuration.satellite = 1;
configuration.carrierToNoiseDensityRatio = 50;
configuration.totalChips = 1023;
configuration.carrierFrequency = 1.57542e9;
configuration.samplingFrequency = 4*1.023e6;
configuration.chippingFrequency = 1.023e6;
configuration.addNoise = true;

phi0             = -2 * pi * configuration.carrierFrequency * (30e3 / 3e8) ; % -2 * pi * fc * (30km / 3e8) I assumed here that 30km is the distance of a satellite to a receiver.
fd               = 0; % Doppler [Hz]
fdr              = 0; % Doppler rate [Hz/s]
doppler_profile  = [phi0, fd, fdr];

configuration.dopplerProfile = doppler_profile;
clear doppler_profile phi0 fd fdr;

save("config_no_doppler.mat", "configuration");

phi0             = -2 * pi * configuration.carrierFrequency * (30e3 / 3e8) ; % -2 * pi * fc * (30km / 3e8) I assumed here that 30km is the distance of a satellite to a receiver.
fd               = 1; % Doppler [Hz]
fdr              = 0; % Doppler rate [Hz/s]
doppler_profile  = [phi0, fd, fdr];

configuration.dopplerProfile = doppler_profile;
clear doppler_profile phi0 fd fdr;

save("config_cte_doppler.mat", "configuration");