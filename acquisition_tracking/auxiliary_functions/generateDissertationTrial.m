function trial = generateDissertationTrial( ...
    baseConfiguration, parameters, seed, carrierToNoiseDensityRatio)
%GENERATEDISSERTATIONTRIAL Build one reproducible causal-channel trial.

rng(seed);
configuration = baseConfiguration;
configuration.carrierToNoiseDensityRatio = carrierToNoiseDensityRatio;

trueDelayEpochRecord = parameters.trueDelay + ...
    cumsum([0 parameters.delayProcessNoiseStd * ...
    randn(1, parameters.simulationSteps)]);
configuration.codeDelay = ...
    repelem(trueDelayEpochRecord, parameters.samplesTotal).';

numberOfCausalTruthTaps = parameters.tapChannelCount;
diffuseTapOrder = 1:(numberOfCausalTruthTaps - 1);
diffusePowerProfile = exp(-0.7 * diffuseTapOrder);
configuration.tdl_channel = zeros(1, numberOfCausalTruthTaps);
configuration.tdl_channel(1) = 1;
configuration.tdl_channel(2:end) = ...
    0.8 * sqrt(diffusePowerProfile) .* ...
    (randn(1, numel(diffuseTapOrder)) + ...
    1j * randn(1, numel(diffuseTapOrder))) / sqrt(2);

[simulatedSignal, ~, ~, LOSDelay] = ...
    gnssReceivedSignal(configuration, parameters.simulationSteps + 1);

truthSampleIndex = round((1:parameters.simulationSteps) * ...
    parameters.samplesTotal);

trial.configuration = configuration;
trial.simulatedSignal = simulatedSignal;
trial.LOSDelay = LOSDelay;
trial.trueDelayRecord = LOSDelay(truthSampleIndex).';
trial.trueDelayEpochRecord = trueDelayEpochRecord;
trial.truthChannelState = configuration.tdl_channel(:);

end
