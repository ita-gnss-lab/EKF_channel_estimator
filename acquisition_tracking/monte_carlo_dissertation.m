clearvars; clc; close all;

scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(genpath(projectRoot));
load(fullfile(scriptDirectory, 'config_no_doppler.mat'));

%% Monte Carlo Setup
baseSeed = 26437226;
numMonteCarloRuns = 2;
numberOfSweepPoints = 7;
monteCarloSeeds = baseSeed + (0:numMonteCarloRuns - 1);

carrierToNoiseSweepDbHz = linspace(20, 55, numberOfSweepPoints);

diagnosticCaseNames = [ ...
    "Adopted"; ...
    "Min delay variance."; ...
    "Max delay variance"; ...
    "Min channel variance"; ...
    "Max channel variance"];

%% Adopted Model Parameters
parameters.simulationSteps = 1000;
parameters.samplesPerChip = 8;
parameters.initialDelayErrorSamples = 0;
parameters.trueDelay = 1e-4;
parameters.nominalCarrierToNoiseDensityRatio = 54.7;
parameters.delayProcessNoiseStd = ...
    0.005 / configuration.chippingFrequency;
parameters.tapProcessNoiseVariance = 5e-3^2;
parameters.initialDelayStd = 0.5 / configuration.chippingFrequency;
parameters.initialChannelStd = 0.5;
parameters.channelSummarySafeguardMs = 500;
parameters.lineWidth = 2;
parameters.fontSize = 13;

parameters.numberOfCorrelators = 17;
parameters.q = (parameters.numberOfCorrelators - 1) / 2;
parameters.C = 2 * parameters.q + 1;
parameters.middleSample = parameters.q + 1;
parameters.tapCount = parameters.q + 1;
parameters.stateDimension = 1 + parameters.tapCount;
parameters.delayStateIndex = 1;
parameters.tapStateIndices = 2:parameters.stateDimension;
parameters.mainTapStateIndex = 2;

configuration.samplingFrequency = ...
    parameters.samplesPerChip * configuration.chippingFrequency;
configuration.applyCarrierPhase = false;
configuration.addNoise = true;
configuration.correlatorHalfSpan = parameters.q;

parameters.epoch = configuration.totalChips / ...
    configuration.chippingFrequency;
parameters.samplesTotal = ...
    round(parameters.epoch * configuration.samplingFrequency);
parameters.channelTapDelays = ...
    (0:parameters.q) / configuration.samplingFrequency;
parameters.timeMs = ...
    (0:parameters.simulationSteps - 1) * parameters.epoch * 1e3;
parameters.chipPeriod = 1 / configuration.chippingFrequency;
parameters.metricStartIndex = ...
    find(parameters.timeMs >= parameters.channelSummarySafeguardMs, 1);
if isempty(parameters.metricStartIndex)
    parameters.metricStartIndex = 1;
end

adoptedDelayProcessNoiseVariance = parameters.delayProcessNoiseStd^2;
adoptedChannelProcessNoiseVariance = parameters.tapProcessNoiseVariance;
delayProcessNoiseVarianceSweep = logspace( ...
    log10((1e-6 / configuration.chippingFrequency)^2), ...
    log10((1e-1 / configuration.chippingFrequency)^2), ...
    numberOfSweepPoints);
channelProcessNoiseVarianceSweep = logspace( ...
    log10(1e-6^2), log10(1e-1^2), numberOfSweepPoints);

totalMonteCarloTasks = ...
    numel(carrierToNoiseSweepDbHz) * numMonteCarloRuns + ...
    2 * numberOfSweepPoints * numMonteCarloRuns + ...
    numMonteCarloRuns + ...
    numel(diagnosticCaseNames);
completedMonteCarloTasks = 0;
monteCarloTimer = tic;

%% Carrier-to-Noise Ratio Sweep
carrierDelayMape = zeros(numel(carrierToNoiseSweepDbHz), numMonteCarloRuns);
carrierChannelMape = zeros(numel(carrierToNoiseSweepDbHz), numMonteCarloRuns);

for carrierIndex = 1:numel(carrierToNoiseSweepDbHz)
    for runIndex = 1:numMonteCarloRuns
        taskName = sprintf("C/N0 sweep %d/%d, run %d/%d", ...
            carrierIndex, numel(carrierToNoiseSweepDbHz), ...
            runIndex, numMonteCarloRuns);
        fprintf("%s\n", taskName);
        trial = generateDissertationTrial( ...
            configuration, parameters, monteCarloSeeds(runIndex), ...
            carrierToNoiseSweepDbHz(carrierIndex));
        result = runDissertationEkf( ...
            trial, parameters, adoptedDelayProcessNoiseVariance, ...
            adoptedChannelProcessNoiseVariance);
        carrierDelayMape(carrierIndex, runIndex) = result.delayMape;
        carrierChannelMape(carrierIndex, runIndex) = result.channelMape;
        completedMonteCarloTasks = printMonteCarloProgress( ...
            taskName, completedMonteCarloTasks, totalMonteCarloTasks, ...
            monteCarloTimer);
    end
end

%% Covariance Sweeps
delayCovarianceDelayMape = zeros(numel(delayProcessNoiseVarianceSweep), ...
    numMonteCarloRuns);
delayCovarianceChannelMape = zeros(numel(delayProcessNoiseVarianceSweep), ...
    numMonteCarloRuns);
channelCovarianceDelayMape = zeros(numel(channelProcessNoiseVarianceSweep), ...
    numMonteCarloRuns);
channelCovarianceChannelMape = zeros(numel(channelProcessNoiseVarianceSweep), ...
    numMonteCarloRuns);
stateCovarianceTraceRecord = zeros(numMonteCarloRuns, ...
    parameters.simulationSteps);

for runIndex = 1:numMonteCarloRuns
    fprintf("Generating adopted-C/N0 trial %d/%d\n", ...
        runIndex, numMonteCarloRuns);
    trial = generateDissertationTrial( ...
        configuration, parameters, monteCarloSeeds(runIndex), ...
        parameters.nominalCarrierToNoiseDensityRatio);

    taskName = sprintf("Adopted state covariance trace, run %d/%d", ...
        runIndex, numMonteCarloRuns);
    fprintf("%s\n", taskName);
    result = runDissertationEkf( ...
        trial, parameters, adoptedDelayProcessNoiseVariance, ...
        adoptedChannelProcessNoiseVariance);
    stateCovarianceTraceRecord(runIndex, :) = ...
        getStateCovarianceTrace(result);
    completedMonteCarloTasks = printMonteCarloProgress( ...
        taskName, completedMonteCarloTasks, totalMonteCarloTasks, ...
        monteCarloTimer);

    for covarianceIndex = 1:numberOfSweepPoints
        taskName = sprintf("Delay covariance sweep %d/%d, run %d/%d", ...
            covarianceIndex, numberOfSweepPoints, ...
            runIndex, numMonteCarloRuns);
        fprintf("%s\n", taskName);
        result = runDissertationEkf( ...
            trial, parameters, ...
            delayProcessNoiseVarianceSweep(covarianceIndex), ...
            adoptedChannelProcessNoiseVariance);
        delayCovarianceDelayMape(covarianceIndex, runIndex) = ...
            result.delayMape;
        delayCovarianceChannelMape(covarianceIndex, runIndex) = ...
            result.channelMape;
        completedMonteCarloTasks = printMonteCarloProgress( ...
            taskName, completedMonteCarloTasks, totalMonteCarloTasks, ...
            monteCarloTimer);

        taskName = sprintf("Channel covariance sweep %d/%d, run %d/%d", ...
            covarianceIndex, numberOfSweepPoints, ...
            runIndex, numMonteCarloRuns);
        fprintf("%s\n", taskName);
        result = runDissertationEkf( ...
            trial, parameters, adoptedDelayProcessNoiseVariance, ...
            channelProcessNoiseVarianceSweep(covarianceIndex));
        channelCovarianceDelayMape(covarianceIndex, runIndex) = ...
            result.delayMape;
        channelCovarianceChannelMape(covarianceIndex, runIndex) = ...
            result.channelMape;
        completedMonteCarloTasks = printMonteCarloProgress( ...
            taskName, completedMonteCarloTasks, totalMonteCarloTasks, ...
            monteCarloTimer);
    end
end

%% Representative Diagnostic Runs
representativeTrial = generateDissertationTrial( ...
    configuration, parameters, monteCarloSeeds(1), ...
    parameters.nominalCarrierToNoiseDensityRatio);

diagnosticDelayCovariances = [ ...
    adoptedDelayProcessNoiseVariance ...
    delayProcessNoiseVarianceSweep(1) ...
    delayProcessNoiseVarianceSweep(end) ...
    adoptedDelayProcessNoiseVariance ...
    adoptedDelayProcessNoiseVariance];
diagnosticChannelCovariances = [ ...
    adoptedChannelProcessNoiseVariance ...
    adoptedChannelProcessNoiseVariance ...
    adoptedChannelProcessNoiseVariance ...
    channelProcessNoiseVarianceSweep(1) ...
    channelProcessNoiseVarianceSweep(end)];
for caseIndex = 1:numel(diagnosticCaseNames)
    taskName = sprintf("Diagnostic %d/%d: %s", ...
        caseIndex, numel(diagnosticCaseNames), ...
        char(diagnosticCaseNames(caseIndex)));
    fprintf("%s\n", taskName);
    diagnosticResults(caseIndex) = runDissertationEkf( ...
        representativeTrial, parameters, ...
        diagnosticDelayCovariances(caseIndex), ...
        diagnosticChannelCovariances(caseIndex));
    completedMonteCarloTasks = printMonteCarloProgress( ...
        taskName, completedMonteCarloTasks, totalMonteCarloTasks, ...
        monteCarloTimer);
end

%% Sweep Plots
plotMeanMapeSweep( ...
    carrierToNoiseSweepDbHz, carrierDelayMape, carrierChannelMape, ...
    "Carrier-to-Noise Ratio Sensitivity", "C/N0 [dB-Hz]", ...
    false, false, parameters.nominalCarrierToNoiseDensityRatio, ...
    parameters);

plotMeanMapeSweep( ...
    delayProcessNoiseVarianceSweep, delayCovarianceDelayMape, ...
    delayCovarianceChannelMape, "Delay Covariance Sensitivity", ...
    "Delay process covariance", true, true, ...
    adoptedDelayProcessNoiseVariance, parameters);

plotMeanMapeSweep( ...
    channelProcessNoiseVarianceSweep, channelCovarianceDelayMape, ...
    channelCovarianceChannelMape, "Channel Covariance Sensitivity", ...
    "Channel process covariance", true, true, ...
    adoptedChannelProcessNoiseVariance, parameters);

%% Diagnostic Plots
plotStateCovarianceTrace(stateCovarianceTraceRecord, parameters);
plotDiagnosticDelayCases(diagnosticResults, diagnosticCaseNames, ...
    parameters);
plotDiagnosticChannelBars(diagnosticResults, diagnosticCaseNames, ...
    parameters);
fprintf("Monte Carlo script finished | elapsed %s\n", ...
    formatElapsedSeconds(toc(monteCarloTimer)));

%% Local Functions
function completedTasks = printMonteCarloProgress( ...
    taskName, completedTasks, totalTasks, timerHandle)

completedTasks = completedTasks + 1;
elapsedSeconds = toc(timerHandle);
remainingSeconds = ...
    elapsedSeconds * (totalTasks - completedTasks) / completedTasks;
estimatedArrival = datestr( ...
    now + remainingSeconds / 86400, 'yyyy-mm-dd HH:MM:SS');

fprintf("%s complete | %d/%d | elapsed %s | ETA %s\n", ...
    char(taskName), completedTasks, totalTasks, ...
    formatElapsedSeconds(elapsedSeconds), estimatedArrival);

end

function timeText = formatElapsedSeconds(elapsedSeconds)

elapsedSeconds = max(0, floor(elapsedSeconds));
hours = floor(elapsedSeconds / 3600);
minutes = floor(mod(elapsedSeconds, 3600) / 60);
secondsValue = mod(elapsedSeconds, 60);
timeText = sprintf("%02d:%02d:%02d", ...
    hours, minutes, secondsValue);

end

function plotMeanMapeSweep( ...
    sweepValues, delayMape, channelMape, figureName, xAxisLabel, ...
    useLogScale, useLogMapeScale, adoptedValue, parameters)

delayMean = mean(delayMape, 2);
delayStd = std(delayMape, 0, 2);
channelMean = mean(channelMape, 2);
channelStd = std(channelMape, 0, 2);
delayColor = [0 0.4470 0.7410];
channelColor = [0.8500 0.1000 0.1000];

figure(Name=figureName, NumberTitle="off");
yyaxis left;
delayHandle = errorbar( ...
    sweepValues, delayMean, delayStd, 'o-', ...
    'LineWidth', parameters.lineWidth, ...
    'Color', delayColor, ...
    'MarkerFaceColor', delayColor);
ylabel("Mean delay MAPE [%]");
ax = gca;
ax.YColor = delayColor;
if useLogScale
    set(ax, "XScale", "log");
end
if useLogMapeScale
    set(ax, "YScale", "log");
end

yyaxis right;
channelHandle = errorbar( ...
    sweepValues, channelMean, channelStd, 's-', ...
    'LineWidth', parameters.lineWidth, ...
    'Color', channelColor, ...
    'MarkerFaceColor', channelColor);
ylabel("Mean channel MAPE [%]");
ax = gca;
ax.YColor = channelColor;
if useLogScale
    set(ax, "XScale", "log");
end
if useLogMapeScale
    set(ax, "YScale", "log");
end
xline(adoptedValue, '--', 'LineWidth', parameters.lineWidth);
setSweepXAxisLimits(sweepValues, useLogScale);
grid on;
title(figureName);
xlabel(xAxisLabel);
set(gca, "FontSize", parameters.fontSize);
legend([delayHandle channelHandle], ...
    {"Delay MAPE", "Channel MAPE"}, "Location", "best");

end

function setSweepXAxisLimits(sweepValues, useLogScale)

xMinimum = min(sweepValues);
xMaximum = max(sweepValues);
if useLogScale
    paddingRatio = (xMaximum / xMinimum)^0.03;
    xlim([xMinimum / paddingRatio xMaximum * paddingRatio]);
else
    padding = 0.03 * (xMaximum - xMinimum);
    if padding == 0
        padding = 1;
    end
    xlim([xMinimum - padding xMaximum + padding]);
end

end

function covarianceTrace = getStateCovarianceTrace(result)

covarianceTrace = zeros(1, numel(result.timeMs));
for k = 1:numel(result.timeMs)
    covarianceTrace(k) = real(trace(result.posteriorCovarianceRecord(:, :, k)));
end

end

function plotStateCovarianceTrace(stateCovarianceTraceRecord, parameters)

meanCovarianceTrace = mean(stateCovarianceTraceRecord, 1);

figure(Name="State Error Covariance", NumberTitle="off");
semilogy(parameters.timeMs, max(meanCovarianceTrace, eps), ...
    'LineWidth', parameters.lineWidth);
grid on;
title("Mean Trace of State Error Covariance");
ylabel("Mean trace of state error covariance");
xlabel("Time [ms]");
set(gca, "FontSize", parameters.fontSize);

end

function plotDiagnosticDelayCases(results, caseNames, parameters)

figure(Name="Diagnostic Delay Estimates", NumberTitle="off");
hold on;
for caseIndex = 1:numel(results)
    plot(parameters.timeMs, ...
        results(caseIndex).delayEstimateRecord / parameters.chipPeriod, ...
        'LineWidth', parameters.lineWidth, ...
        'DisplayName', caseNames(caseIndex));
end
plot(parameters.timeMs, ...
    results(1).trueDelayRecord / parameters.chipPeriod, ...
    '--k', 'LineWidth', 1.5 * parameters.lineWidth, ...
    'DisplayName', "Truth");
hold off;
grid on;
title("Diagnostic Delay Estimates");
ylabel("LOS Delay [T_c]");
xlabel("Time [ms]");
legend("Location", "best");
set(gca, "FontSize", parameters.fontSize);

end

function plotDiagnosticChannelBars(results, caseNames, parameters)

channelTapNumbers = 0:(parameters.tapCount - 1);
numberOfCases = numel(results);
estimatedMean = zeros(parameters.tapCount, numberOfCases);
estimatedRealStd = zeros(parameters.tapCount, numberOfCases);
estimatedImagStd = zeros(parameters.tapCount, numberOfCases);

for caseIndex = 1:numberOfCases
    metricSamples = parameters.metricStartIndex:parameters.simulationSteps;
    channelSummaryRecord = ...
        results(caseIndex).channelStateRecord(:, metricSamples);
    estimatedMean(:, caseIndex) = mean(channelSummaryRecord, 2);
    estimatedRealStd(:, caseIndex) = std(real(channelSummaryRecord), 0, 2);
    estimatedImagStd(:, caseIndex) = std(imag(channelSummaryRecord), 0, 2);
end

trueChannel = results(1).truthChannelState;

figure(Name="Diagnostic Channel Coefficient Bars", NumberTitle="off");
tiledlayout(1, 2);

plotDiagnosticCoefficientBar( ...
    channelTapNumbers, real(estimatedMean), real(trueChannel), ...
    estimatedRealStd, "Real Channel Coefficients", caseNames, parameters);
plotDiagnosticCoefficientBar( ...
    channelTapNumbers, imag(estimatedMean), imag(trueChannel), ...
    estimatedImagStd, "Imaginary Channel Coefficients", caseNames, ...
    parameters);

end

function plotDiagnosticCoefficientBar( ...
    tapNumbers, estimatedValues, trueValues, errorValues, plotTitle, ...
    caseNames, parameters)

nexttile;
numberOfCases = size(estimatedValues, 2);
numberOfSeries = numberOfCases + 1;
tapPositions = 2.1 * (1:numel(tapNumbers));
barOffsets = linspace(-0.62, 0.62, numberOfSeries);
barWidth = 0.14;
estimateColors = lines(numberOfCases);
barHandles = gobjects(numberOfSeries, 1);

hold on;
for caseIndex = 1:numberOfCases
    xPositions = tapPositions + barOffsets(caseIndex);
    barHandles(caseIndex) = bar( ...
        xPositions, estimatedValues(:, caseIndex), barWidth);
    barHandles(caseIndex).FaceColor = estimateColors(caseIndex, :);
    errorbar( ...
        xPositions, estimatedValues(:, caseIndex), ...
        errorValues(:, caseIndex), ...
        'LineStyle', 'none', ...
        'Color', estimateColors(caseIndex, :), ...
        'LineWidth', parameters.lineWidth);
end

truthPositions = tapPositions + barOffsets(end);
barHandles(end) = bar(truthPositions, trueValues, barWidth);
barHandles(end).FaceColor = [0 0.5 0];
hold off;
grid on;
legend(barHandles, [caseNames(:); "Truth"], "Location", "eastoutside");
title(plotTitle);
ylabel("Channel Coefficients");
xticks(tapPositions);
xticklabels(compose("h_{%d}", tapNumbers));
xlim([tapPositions(1) - 1.1 tapPositions(end) + 1.1]);
set(gca, "FontSize", parameters.fontSize);

end
