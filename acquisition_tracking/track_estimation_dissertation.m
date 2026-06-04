clearvars; clc; close all;

scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(genpath(projectRoot));
load(fullfile(scriptDirectory, 'config_no_doppler.mat'));

%% Simulation Setup
rngSeed = 26437226;
simulationSteps = 2000;

% Andreas Iliopoulos article, Section 4: lambda = 4 and
% fs = 2*lambda/Tc, which gives 8 samples per chip.
samplesPerChip = 8;

% Andreas Iliopoulos article, Section 4: initial delay estimation error is
% tau_0^(0) - tau_0^(0)+ = 0.0*Tc.
initialDelayErrorSamples = 0;

% (Rodrigo): I configured the absolute simulation delay reference as a satellite away by 30km.
trueDelay = 1e-4;

configuration.samplingFrequency = ...
    samplesPerChip * configuration.chippingFrequency;

% Andreas Iliopoulos article, Section 4: nominal C/N0 = 54.7 dB-Hz.
configuration.carrierToNoiseDensityRatio = 54.7;

% (Rodrigo): I configured the simplified no-Doppler/no-carrier-phase
% scenario.
configuration.applyCarrierPhase = false;
configuration.addNoise = true;

%% Multicorrelator Setup
% Andreas Iliopoulos article, Section 4: P = 2L + 1 = 17 correlators.
numberOfCorrelators = 17;
correlatorHalfSpan = (numberOfCorrelators - 1) / 2;
q = correlatorHalfSpan;
C = 2 * q + 1;
middleSample = q + 1;
tapCount = q + 1;

configuration.correlatorHalfSpan = q;

epoch = configuration.totalChips / configuration.chippingFrequency;
samplesTotal = round(epoch * configuration.samplingFrequency);

% Andreas Iliopoulos article, Section 4: correlator spacing
% Delta = 2*Tc/L = 0.125*Tc. With fs = 8/Tc, one sample is 0.125*Tc.
channelTapDelays = (0:q) / configuration.samplingFrequency;

%% EKF Model Parameters
stateDimension = 1 + tapCount;
delayStateIndex = 1;
tapStateIndices = 2:stateDimension;
mainTapStateIndex = 2;

% (Rodrigo): I configured the simulated delay random walk. Andreas uses
% beta_tau = 0.001*Tc in Section 4; this script uses 0.005*Tc. The system
% only works robustly with this configuration below. With Andreas'
% configuration, it only works for some time (he used 200 ms/epochs in his
% work), but later it loses lock.
delayProcessNoiseStd = 0.005 / configuration.chippingFrequency;

% Andreas Iliopoulos article, Section 4: initial LOS delay standard
% deviation is 0.5*Tc.
initialDelayStd = 0.5 / configuration.chippingFrequency;

% (Rodrigo): I configured the channel process noise for this dissertation
% experiment. Andreas uses beta_h = 0.01 for the channel random walk in
% Section 4. Larger values can degrade the delay estimation precision.
tapProcessNoiseVariance = 5e-3^2;

% (Rodrigo): I configured this initial channel covariance. Andreas states an
% initial amplitude standard deviation of 0.1 in Section 4. I configured
% with as 0.5, because I noticed that the channel taps converges faster to
% the true values this way.
initialChannelStd = 0.5;

% (Rodrigo): I configured a 500 ms safeguard before computing the channel
% coefficient mean and standard deviation, so the initial convergence
% transient is not included in the summary bars.
channelSummarySafeguardMs = 500;

epochVector = 1:simulationSteps;
timeMs = (epochVector - 1) * epoch * 1e3;
chipPeriod = 1 / configuration.chippingFrequency;
channelSummaryStartIndex = find(timeMs >= channelSummarySafeguardMs, 1);
if isempty(channelSummaryStartIndex)
    channelSummaryStartIndex = 1;
end

parameters.simulationSteps = simulationSteps;
parameters.samplesPerChip = samplesPerChip;
parameters.initialDelayErrorSamples = initialDelayErrorSamples;
parameters.trueDelay = trueDelay;
parameters.delayProcessNoiseStd = delayProcessNoiseStd;
parameters.tapProcessNoiseVariance = tapProcessNoiseVariance;
parameters.initialDelayStd = initialDelayStd;
parameters.initialChannelStd = initialChannelStd;
parameters.numberOfCorrelators = numberOfCorrelators;
parameters.q = q;
parameters.C = C;
parameters.middleSample = middleSample;
parameters.tapCount = tapCount;
parameters.stateDimension = stateDimension;
parameters.delayStateIndex = delayStateIndex;
parameters.tapStateIndices = tapStateIndices;
parameters.mainTapStateIndex = mainTapStateIndex;
parameters.epoch = epoch;
parameters.samplesTotal = samplesTotal;
parameters.channelTapDelays = channelTapDelays;
parameters.timeMs = timeMs;
parameters.chipPeriod = chipPeriod;
parameters.metricStartIndex = channelSummaryStartIndex;

%% Synthetic Channel, Delay Truth, and EKF Simulation
% (Rodrigo): I configured this synthetic causal diffuse channel model in
% generateDissertationTrial. It is not the analog distortion-fault model
% from Figure 12 of the Andreas article. The idea behind this is to model
% the specular component (LOS) as 1 and a diffuse complex channel profile
% with an exponential decaying structure. Note that this is purely
% synthetic and would not necessarily behave as a GNSS-R measurement
% channel.
trial = generateDissertationTrial( ...
    configuration, parameters, rngSeed, configuration.carrierToNoiseDensityRatio);
result = runDissertationEkf( ...
    trial, parameters, delayProcessNoiseStd^2, tapProcessNoiseVariance);

configuration = trial.configuration;
trueDelayEpochRecord = trial.trueDelayEpochRecord;
truthChannelState = result.truthChannelState;
trueDelayRecord = result.trueDelayRecord;
delayEstimateRecord = result.delayEstimateRecord;
channelStateRecord = result.channelStateRecord;
innovationRecord = result.innovationRecord;

%% Plot Preparation
% (Rodrigo): I configured the plotting style; these values are not
% simulation parameters from the Andreas Iliopoulos article.
lineWidth = 2;
fontSize = 13;

numberOfChannelTaps = tapCount;
channelTapNumbers = 0:(numberOfChannelTaps - 1);
channelTapDelayTc = channelTapNumbers / samplesPerChip;
channelPlotRows = ceil(sqrt(numberOfChannelTaps));
channelPlotColumns = ceil(numberOfChannelTaps / channelPlotRows);

channelSummaryRecord = channelStateRecord(:, channelSummaryStartIndex:end);

%% Delay and Innovation Figures
figure(Name="STD of the innovations", NumberTitle="off");
plot(timeMs, std(innovationRecord, 1, 1), 'LineWidth', lineWidth);
hold on;
plot(timeMs, zeros(1, simulationSteps), '--', 'LineWidth', lineWidth);
hold off;
legend({"Innovation STD", "Truth"});
ylabel("Standard deviation of the innovations");
xlabel("Time [ms]");
set(gca, "FontSize", fontSize);

figure(Name="Middle tap of the innovation sequence", NumberTitle="off");
plot(timeMs, real(innovationRecord(middleSample, :)), ...
    'LineWidth', lineWidth);
hold on;
plot(timeMs, imag(innovationRecord(middleSample, :)), ...
    'LineWidth', lineWidth);
plot(timeMs, zeros(1, simulationSteps), '--', 'LineWidth', lineWidth);
hold off;
legend({"Real", "Imaginary", "Truth"});
ylabel("Innovation sequence of the middle tap");
xlabel("Time [ms]");
set(gca, "FontSize", fontSize);

figure(Name="Delay Estimation", NumberTitle="off");
plot(timeMs, delayEstimateRecord / chipPeriod, 'LineWidth', lineWidth);
hold on;
plot(timeMs, trueDelayRecord / chipPeriod, '--', 'LineWidth', lineWidth);
hold off;
legend({"EKF's estimated delay", "True delay"});
ylabel("LOS Delay [T_c]");
xlabel("Time [ms]");
set(gca, "FontSize", fontSize);

%% Channel Time History Figures
plotTapHistories( ...
    "Real Channel Tap Estimates", timeMs, ...
    real(channelStateRecord), real(truthChannelState), ...
    channelTapNumbers, channelPlotRows, channelPlotColumns, ...
    "Real", lineWidth, fontSize);

plotTapHistories( ...
    "Imaginary Channel Tap Estimates", timeMs, ...
    imag(channelStateRecord), imag(truthChannelState), ...
    channelTapNumbers, channelPlotRows, channelPlotColumns, ...
    "Imaginary", lineWidth, fontSize);

%% Channel Summary Figures
plotCoefficientSummary( ...
    channelTapNumbers, channelSummaryRecord, truthChannelState, ...
    lineWidth, fontSize);

plotImpulseResponseSurfaces( ...
    timeMs, channelTapDelayTc, channelStateRecord, ...
    truthChannelState, simulationSteps, fontSize);

plotComplexTapScatter( ...
    channelSummaryRecord, truthChannelState, channelTapNumbers, ...
    lineWidth, fontSize);

%% Plot Functions
function plotTapHistories( ...
    figureName, timeMs, estimatedTaps, trueTaps, tapNumbers, ...
    plotRows, plotColumns, yAxisLabel, lineWidth, fontSize)

figure(Name=figureName, NumberTitle="off");
tiledlayout(plotRows, plotColumns);
for tapIndex = 1:numel(tapNumbers)
    nexttile;
    plot(timeMs, estimatedTaps(tapIndex, :), 'LineWidth', lineWidth);
    hold on;
    plot(timeMs, trueTaps(tapIndex) * ones(size(timeMs)), ...
        '--', 'LineWidth', lineWidth);
    hold off;
    title(sprintf("Tap %+d", tapNumbers(tapIndex)));
    ylabel(yAxisLabel);
    xlabel("Time [ms]");
    set(gca, "FontSize", fontSize);
end
legend({"Estimate", "Truth"});

end

function plotCoefficientSummary( ...
    tapNumbers, estimatedChannel, trueChannel, lineWidth, fontSize)

estimatedMean = mean(estimatedChannel, 2);
estimatedRealStd = std(real(estimatedChannel), 0, 2);
estimatedImagStd = std(imag(estimatedChannel), 0, 2);

figure(Name="Channel Coefficient Summary", NumberTitle="off");
tiledlayout(1, 2);

plotCoefficientBar( ...
    tapNumbers, real(estimatedMean), real(trueChannel), estimatedRealStd, ...
    "Real Channel Coefficients", lineWidth, fontSize);
plotCoefficientBar( ...
    tapNumbers, imag(estimatedMean), imag(trueChannel), estimatedImagStd, ...
    "Imaginary Channel Coefficients", lineWidth, fontSize);

end

function plotCoefficientBar( ...
    tapNumbers, estimatedValues, trueValues, errorValues, ...
    plotTitle, lineWidth, fontSize)

nexttile;
tapPositions = 1.5 * (1:numel(tapNumbers));
estimatedPositions = tapPositions - 0.22;
truthPositions = tapPositions + 0.22;

estimatedBar = bar(estimatedPositions, estimatedValues, 0.32);
estimatedBar.FaceColor = [0.55 0.75 0.95];
hold on;
truthBar = bar(truthPositions, trueValues, 0.32);
truthBar.FaceColor = [0 0.5 0];
errorbar(estimatedPositions, estimatedValues, errorValues, ...
    'LineStyle', 'none', 'Color', [0 0.4470 0.7410], ...
    'LineWidth', lineWidth);
hold off;
grid on;
legend([estimatedBar truthBar], ...
    {"Estimated Channel Coefficients", "True Channel Coefficients"}, ...
    "Location", "best");
title(plotTitle);
ylabel("Channel Coefficients");
xticks(tapPositions);
xticklabels(compose("h_{%d}", tapNumbers));
xlim([tapPositions(1) - 0.8 tapPositions(end) + 0.8]);
set(gca, "FontSize", fontSize);

end

function plotImpulseResponseSurfaces( ...
    timeMs, tapDelayTc, estimatedChannel, trueChannel, ...
    simulationSteps, fontSize)

[timeGrid, tapDelayGrid] = meshgrid(timeMs, tapDelayTc);

figure(Name="Channel Impulse Response History", NumberTitle="off");
tiledlayout(1, 2);

plotImpulseSurfaceTile( ...
    timeGrid, tapDelayGrid, real(estimatedChannel), ...
    repmat(real(trueChannel), 1, simulationSteps), ...
    "Real Channel Impulse Response", "Real", fontSize);
plotImpulseSurfaceTile( ...
    timeGrid, tapDelayGrid, imag(estimatedChannel), ...
    repmat(imag(trueChannel), 1, simulationSteps), ...
    "Imaginary Channel Impulse Response", "Imaginary", fontSize);

end

function plotImpulseSurfaceTile( ...
    timeGrid, tapDelayGrid, estimatedSurface, trueSurface, ...
    plotTitle, zAxisLabel, fontSize)

nexttile;
surf(timeGrid, tapDelayGrid, estimatedSurface, 'EdgeColor', 'none');
hold on;
surf(timeGrid, tapDelayGrid, trueSurface, ...
    'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.85, ...
    'EdgeColor', 'none');
hold off;
grid on;
colormap(gca, turbo);
colorbar;
view(42, 28);
camproj("perspective");
axis tight;
legend({"Estimate", "Truth"});
title(plotTitle);
xlabel("Time [ms]");
ylabel("Delay [T_c]");
zlabel(zAxisLabel);
set(gca, "FontSize", fontSize);

end

function plotComplexTapScatter( ...
    estimatedChannel, trueChannel, tapNumbers, lineWidth, fontSize)

figure(Name="Final Channel Tap Estimates", NumberTitle="off");
plot(real(trueChannel), imag(trueChannel), 'x', ...
    'LineWidth', lineWidth, 'MarkerSize', 10);
hold on;
plot(real(estimatedChannel(:, end)), imag(estimatedChannel(:, end)), ...
    'o', 'LineWidth', lineWidth, 'MarkerSize', 8);
for tapIndex = 1:numel(tapNumbers)
    text(mean(real(estimatedChannel(tapIndex, :))), ...
        mean(imag(estimatedChannel(tapIndex, :))), ...
        sprintf(" %+d", tapNumbers(tapIndex)));
end
hold off;
grid on;
axis equal;
legend({"Truth", "Final estimate"});
xlabel("Real");
ylabel("Imaginary");
set(gca, "FontSize", fontSize);

end
