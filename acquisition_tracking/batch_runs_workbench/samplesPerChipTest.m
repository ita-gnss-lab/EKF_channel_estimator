%% Default simulation variables 

simulation_config.satellite =                   1;
simulation_config.carrierToNoiseDensityRatio =  50; %dB
simulation_config.totalChips =                  1023;
simulation_config.carrierFrequency =            1.57542e+09; %Hz
simulation_config.chippingFrequency =           1023000; %Hz
simulation_config.samplesPerChip =              16;
simulation_config.addNoise =                    true; 
simulation_config.dopplerProfile = ...
[-2 * pi * simulation_config.carrierFrequency * (30e3 / 3e8) ...
    0 0];
simulation_config.correlatorHalfSpaan =         5;


%% Default filter variables

filter_config.numberOfCausalTruthTaps    = 9;
filter_config.initialDelayErrorSamples   = 0;
filter_config.phaseErrorStd              = 0;

%% Default seed

seed = 26437226;

%% Batch run definition
% Set the experiment prefix
experimentPrefix = 'samplesPerChipTest';

% Preset the batch run
batchVector = 16:8:64;

% In the loop, preset the cahnging variable before run_simulation
for i = 1:length(batchVector)
    simulation_config.samplesPerChip = batchVector(i);
    filter_config.initialDelayErrorSamples   = 0;
    experimentName = [experimentPrefix  '_'  int2str(batchVector(i)) 'dB'];
    fprintf('%s running\n Progress: Experiment %d/%d\n', experimentName, i, length(batchVector));
    run_simulation(experimentName, simulation_config, filter_config, seed);
end
