load config_cte_doppler.mat
configuration.addNoise = true; 
noiseVector = [50 45 40 35 30 25];
nominalDoppler = 2*pi*configuration.dopplerProfile(2);
threshold = 2*pi*300;
myExitCondition = @(state) abs(state(3) - nominalDoppler) > threshold;
saveTo = 'noiseRobustnessTest';

initialCondition.mainTap = 1;
initialCondition.phaseError = 0;
initialCondition.DopplerError = 50;
initialCondition.delay = 1.005e-4;

for i = 1:length(noiseVector)
    configuration.carrierToNoiseDensityRatio = noiseVector(i);
    testbench_function(configuration, i, myExitCondition, saveTo, initialCondition);
end




