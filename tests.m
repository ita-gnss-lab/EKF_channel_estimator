clearvars; clc;

% Set the path to all files inside this repository.
addpath(genpath(fullfile("..","EKF_channel_estimator")));

load(fullfile("gnss_simulator","config_no_doppler.mat"))
rng(26437226);

%% Parameters
simulationSteps = 500;
numberOfTaps = 2;
% NOTE: Epoch is equivalent to T = (N_c / (1/T_c)) in our work.
epoch = configuration.totalChips / configuration.chippingFrequency;
% NOTE: What this would represent?
carrierError = 1:4;
% NOTE: This is outputing a 1x0 empty double row vector.
channelWeights = [5, 0:(numberOfTaps - 1)];


%% Simulate Signal
[simulatedSignal, totalTime] = gnss_received_signal(configuration, simulationSteps + 1);
% NOTE: samplesTotal = N in our notation, i.e., the amount of samples 
% within a PR-code block.
samplesTotal = epoch*configuration.samplingFrequency;

%% Checking if the simulated signal can be acquired
% NOTE: It seems that the signal is being acquired correctly.
gsa = gnssSignalAcquirer( ...
    "SampleRate", configuration.samplingFrequency, ...
    "GNSSSignalType", "GPS C/A", ...
    "FrequencyRange", [-10e3, 10e3], ...
    "FrequencyResolution", 500);

% NOTE: If you change the line below to simulationSteps/2, it becomes clear
% that the code-phase offset is changing over time, as you also have shown
% in your section below, where you plot the correlation of the first
% PR-code block with shifted samples over time.
samplesOffset = 1000;
[acqtable, corrmat] = gsa(simulatedSignal(samplesOffset + 1: (samplesOffset + configuration.samplingFrequency*1e-3)), 1);
disp("Code-Phase Offset: " + num2str(acqtable.CodePhaseOffset));
disp("Coarse Doppler Shift: " + num2str(acqtable.FrequencyOffset));
freqRange= gsa.FrequencyRange; % Range of the frequency search in Hz
stepSize = gsa.FrequencyResolution; % Step size of frequency search in Hz
satIndex = 1; 
figure;
% Visualize 1st satellite correlation
mesh(freqRange(1):stepSize:freqRange(2),0:size(corrmat,1)-1,corrmat(:,:,satIndex)) % Surface plot
xlabel("Doppler Offset (Hz)");
ylabel("Code Phase Offset (samples)");
zlabel("Correlation");
title("Correlation Plot for PRN ID: " + acqtable.PRNID(satIndex));

%%
CorrelationResults = zeros(simulationSteps, 1);
for i = 1:(simulationSteps - 2)
    CorrelationResults(i) = abs(simulatedSignal((samplesTotal + 1):(2)*samplesTotal)).' * ...
    abs(simulatedSignal(((i+1)*samplesTotal + 1):(i+2)*samplesTotal)) / samplesTotal;
end

figure;
plot(CorrelationResults);

in = (1:4000)';
delay_handler = dsp.VariableFractionalDelay("InterpolationMethod","Linear", 'MaximumDelay',9999);
delayVec = 100*ones(4000, 1);
outcase1 = delay_handler(in,delayVec);
start = circshift(in, delayVec);
%outcase1(1:delayVec) = start(1:delayVec);
outcase1(delayVec) = start(delayVec);

figure;
plot(outcase1(:,1));
hold on;
plot(in)
hold off;