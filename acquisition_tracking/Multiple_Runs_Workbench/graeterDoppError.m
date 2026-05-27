load config_cte_doppler.mat
configuration.addNoise = false; 
nominalDoppler = 2*pi*configuration.dopplerProfile(2);
threshold = 2*pi*300;
myExitCondition = @(state) abs(state(3) - nominalDoppler) > threshold;
saveTo = 'mainTapRobustnessTest';

initialCondition.mainTap = 1;
initialCondition.phaseError = 0;
initialCondition.DopplerError = 50;
initialCondition.delay = 1.005e-4;
doppErrorVector = 5:15:50

for i = 1:length(tapVector)
    initialCondition.mainTap = tapVector(i);
    testbench_function(configuration, i, myExitCondition, saveTo, initialCondition);
end