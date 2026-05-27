load config_cte_doppler.mat
configuration.addNoise = false; 
nominalDoppler = 2*pi*configuration.dopplerProfile(2);
threshold = 0.0101;
myExitCondition = @(state) abs(state(4+5+1) - 1) > threshold;
saveTo = 'teste_simple_system_4';

initialCondition.mainTap = 1.01;
initialCondition.phaseError = pi/4;
initialCondition.DopplerError = 0;
initialCondition.delay = 1.005e-4;

varr = -10: 1 : -5;
varr = 10 .^varr;
results = zeros(1, length(varr));

for i = 1:length(varr)
    results(1, i) = testbench_function_simple(configuration, i, myExitCondition, saveTo, initialCondition, varr(i));
end

