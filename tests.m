clear;
load config_no_doppler.mat
rng(26437226);

%% Parameters
simulationSteps = 500;
numberOfTaps = 2;
epoch = configuration.totalChips / configuration.chippingFrequency;
carrierError = 1:4;
channelWeights = 5 + 0:numberOfTaps;


%% Simulate Signal
% (Rodrigo): Put this out of the loop
[simulatedSignal, totalTime] = gnss_received_signal(configuration, simulationSteps + 1);
samplesTotal = epoch*configuration.samplingFrequency + 1;

result = zeros(simulationSteps, 1);
for i = 1:(simulationSteps - 2)
    result(i) = real(simulatedSignal((samplesTotal + 1):(2)*samplesTotal)).' * ...
    real(simulatedSignal(((i+1)*samplesTotal + 1):(i+2)*samplesTotal)) / samplesTotal;
end

plot(result);


in = (1:4000)';
delay_handler = dsp.VariableFractionalDelay("InterpolationMethod","Linear", 'MaximumDelay',9999);
delayVec = 100*ones(4000, 1);
outcase1 = delay_handler(in,delayVec);
start = circshift(in, delayVec);
outcase1(1:delayVec) = start(1:delayVec);

plot(outcase1(:,1));
hold on;
plot(in)
hold off;