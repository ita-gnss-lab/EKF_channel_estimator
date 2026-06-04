clearvars; clc; close all;

addpath(genpath(fullfile("..", "..","EKF_channel_estimator")));

load config_no_doppler.mat
rng(26437226);
samplesPerChip = 8;
configuration.samplingFrequency = samplesPerChip * configuration.chippingFrequency;

%% Simulation Setup
simulationSteps = 500;
configuration.carrierToNoiseDensityRatio = 54.7;
trueDelay = 1e-4;
configuration.applyCarrierPhase = false;

numberOfCorrelators = 17;
correlatorHalfSpan = (numberOfCorrelators - 1) / 2;
numberOfCausalTruthTaps = correlatorHalfSpan + 1;
diffuseTapOrder = 1:(numberOfCausalTruthTaps - 1);
diffusePowerProfile = exp(-0.7 * diffuseTapOrder);
configuration.tdl_channel = zeros(1, numberOfCausalTruthTaps);
configuration.tdl_channel(1) = 1;
configuration.tdl_channel(2:end) = 0.8 * sqrt(diffusePowerProfile) .* ...
    (randn(1, numel(diffuseTapOrder)) + ...
    1j * randn(1, numel(diffuseTapOrder))) / sqrt(2);

%% Configuration Flags
useTapEnergyConstraint = false;
useAcausalTaps = false;
usePerfectFrozenTruthState = false;
isAdaptiveMeasurementCovariance = false;
isFixedDelayJacobian = true;
isAdaptiveStateCovariance = false;
configuration.addNoise = true;
plotMeasures = false;

%% Parameters
constraint_noise = 10^(-2.63);
q = correlatorHalfSpan;
C = 2*q + 1;
middleSample = q + 1;
tapCount = 2*q + 1;
stateDimension = 1 + tapCount;
delayStateIndex = 1;
tapStateIndices = 2:stateDimension;
mainTapStateIndex = q + 2;
acausalTapStateIndices = 2:(q + 1);
acausalTapChannelIndices = 1:q;
epoch = configuration.totalChips / configuration.chippingFrequency;
configuration.correlatorHalfSpan = q;
samplesTotal = round(epoch * configuration.samplingFrequency);
delayProcessNoiseStd = 0.005 / configuration.chippingFrequency;
trueDelayEpochRecord = trueDelay + ...
    cumsum([0 delayProcessNoiseStd * randn(1, simulationSteps)]);
configuration.codeDelay = repelem(trueDelayEpochRecord, samplesTotal).';
channelTapDelays = (-q:q) / configuration.samplingFrequency;
delayJacobianStep = 1e-3 / configuration.samplingFrequency;
adaptiveStateCovarianceWindow = 50;
initialDelayErrorSamples = 0;
initialDelayEstimate = trueDelayEpochRecord(1) + ...
    initialDelayErrorSamples / configuration.samplingFrequency;
initialDelayStd = 0.5 / configuration.chippingFrequency;

truthChannelAcausal = zeros(tapCount, 1);
numberOfTruthTaps = min(numel(configuration.tdl_channel), q + 1);
truthChannelAcausal(q + 1:q + numberOfTruthTaps) = ...
    configuration.tdl_channel(1:numberOfTruthTaps).';
truthChannelState = truthChannelAcausal;

%% Covariances
carrierToNoiseRatioLinear = 10^(configuration.carrierToNoiseDensityRatio / 10);
thermalNoiseVariance = configuration.samplingFrequency / carrierToNoiseRatioLinear;

delayProcessNoiseVariance = delayProcessNoiseStd^2;
tapProcessNoiseVariance = 1e-4^2;
channelProcessCovariance = tapProcessNoiseVariance * eye(tapCount);
channelProcessCovariance(q + 1, q + 1) = ...
    10 * channelProcessCovariance(q + 1, q + 1);
QBase = blkdiag(delayProcessNoiseVariance, channelProcessCovariance);
if ~useAcausalTaps
    QBase(acausalTapStateIndices, :) = 0;
    QBase(:, acausalTapStateIndices) = 0;
end
if usePerfectFrozenTruthState
    QBase(:) = 0;
end
Q = QBase;

%% State History Vectors
delayEstimateRecord = zeros(1, simulationSteps);
channelStateRecord = zeros(tapCount, simulationSteps);
measurementDimension = C + double(useTapEnergyConstraint);
innovationRecord = zeros(measurementDimension, simulationSteps);
kalmanGainRecord = zeros(stateDimension, measurementDimension, simulationSteps);
constraintRecord = zeros(1, simulationSteps);
adaptiveStateCovarianceErrorMemory = ...
    zeros(measurementDimension, adaptiveStateCovarianceWindow);

%% Initialization
x_k_k_1 = zeros(stateDimension, 1);
x_k_k_1(delayStateIndex) = initialDelayEstimate;
x_k_k_1(mainTapStateIndex) = 1;

otherTapMask = true(size(x_k_k_1));
otherTapMask([delayStateIndex mainTapStateIndex]) = false;

if usePerfectFrozenTruthState
    x_k_k_1(delayStateIndex) = trueDelayEpochRecord(1);
    x_k_k_1(tapStateIndices) = truthChannelState;
end

initialChannelCovarianceMatrix = 0.5^2 * eye(tapCount);
if usePerfectFrozenTruthState
    initialChannelCovarianceMatrix(:) = 0;
end
if ~useAcausalTaps
    initialChannelCovarianceMatrix(acausalTapChannelIndices, :) = 0;
    initialChannelCovarianceMatrix(:, acausalTapChannelIndices) = 0;
end
P_k_k_1 = blkdiag(initialDelayStd^2, initialChannelCovarianceMatrix);

correlatorBank = buildCorrelatorBank( ...
    configuration, x_k_k_1(delayStateIndex), q);
R = (thermalNoiseVariance / samplesTotal) * ...
    ((correlatorBank * correlatorBank') / samplesTotal);
if useTapEnergyConstraint
    R = [R zeros(C, 1); zeros(1, C) constraint_noise];
end

%% Simulate Signal
[simulatedSignal, ~, ~, LOSDelay] = gnssReceivedSignal(configuration, simulationSteps + 1);

%% Simulation
correlatorTaps = -q:1:q;
for k = 1:simulationSteps
    receivedSignal = simulatedSignal(((k - 1) * samplesTotal + 1: k * samplesTotal));

    if k > 1
        %% EKF Update Step
        correlatorBank = buildCorrelatorBank( ...
            configuration, x_k_k_1(delayStateIndex), q);
        if isAdaptiveMeasurementCovariance
            R = (thermalNoiseVariance / samplesTotal) * ...
                ((correlatorBank * correlatorBank') / samplesTotal);
            if useTapEnergyConstraint
                R = [R zeros(C, 1); zeros(1, C) constraint_noise];
            end
        end

        z_k = (correlatorBank * receivedSignal) / samplesTotal;
        channelWeights = x_k_k_1(tapStateIndices);
        currentShiftedCorrelations = getShiftedCorrelationsFromBank( ...
            correlatorBank, ...
            x_k_k_1(delayStateIndex), ...
            configuration, ...
            channelTapDelays);
        channelWeightsJacobian = ...
            currentShiftedCorrelations / samplesTotal;
        z_hat_k_aux = channelWeightsJacobian * channelWeights;
        z_hat_k = z_hat_k_aux;
        if useTapEnergyConstraint
            mainTapValue = x_k_k_1(mainTapStateIndex);
            mainTapEnergy = abs(mainTapValue)^2;
            otherTapEnergy = sum(abs(x_k_k_1(otherTapMask)).^2);
            constraint_value = ...
                otherTapEnergy / ((q - 1) * mainTapEnergy);
            z_k = [z_k; 0];
            z_hat_k = [z_hat_k; constraint_value];
            constraintRecord(:, k) = constraint_value;
        end

        if plotMeasures
            plot(correlatorTaps, real(z_k(1:C)));
            hold on;
            plot(correlatorTaps, real(z_hat_k(1:C)));
            plot(correlatorTaps, imag(z_k(1:C)));
            plot(correlatorTaps, imag(z_hat_k(1:C)));
            plot(correlatorTaps(q+1:end), real(z_k(q+1:C)), 'o');
            plot(correlatorTaps(q+1:end), real(z_hat_k(q+1:C)), 'x');
            hold off;
            ylabel('Real and Imag parts of z_k and z_k_hat');
            xlabel(sprintf('Correlator tap %d', k));
            legend({'Real $z[k]$', 'Real $\hat{z}[k]$', ...
                'Imag $z[k]$', 'Imag $\hat{z}[k]$', ...
                'Causal real $z[k]$', 'Causal real $\hat{z}[k]$'}, ...
                'Interpreter','latex');
            pause(0.01);
        end

        if isFixedDelayJacobian
            delayJacobian = delayJacobianFunctionSimplified_acausal( ...
                0, ...
                channelWeights, ...
                1 / configuration.samplingFrequency, ...
                q, ...
                1 / configuration.chippingFrequency, ...
                1);
        else
            delayPlusCorrelations = getShiftedCorrelationsFromBank( ...
                correlatorBank, ...
                x_k_k_1(delayStateIndex) + delayJacobianStep, ...
                configuration, ...
                channelTapDelays);
            delayMinusCorrelations = getShiftedCorrelationsFromBank( ...
                correlatorBank, ...
                x_k_k_1(delayStateIndex) - delayJacobianStep, ...
                configuration, ...
                channelTapDelays);
            delayJacobian = ...
                ((delayPlusCorrelations - delayMinusCorrelations) * ...
                channelWeights) / (2 * delayJacobianStep * samplesTotal);
        end

        jacobian = [delayJacobian channelWeightsJacobian];
        if useTapEnergyConstraint
            losParcel = ...
                -mainTapValue * otherTapEnergy / mainTapEnergy;
            tapsParcel = ...
                (2 / ((q - 1) * mainTapEnergy)) * ...
                [x_k_k_1(2:q+1); losParcel; x_k_k_1(q+3:end)];
            constraintLine = [0 tapsParcel'];
            jacobian = [jacobian; constraintLine];
        end

        PJacobianTranspose = P_k_k_1 * jacobian';
        innovationCovariance = jacobian * PJacobianTranspose + R;
        K_k = PJacobianTranspose * ...
            (innovationCovariance \ eye(measurementDimension));
        kalmanGainRecord(:, :, k) = K_k;

        innovation = z_k - z_hat_k;
        innovationRecord(:, k) = innovation;

        x_k_k = x_k_k_1 + K_k * innovation;
        x_k_k(delayStateIndex) = real(x_k_k(delayStateIndex));
        if ~useAcausalTaps
            x_k_k(acausalTapStateIndices) = 0;
        end

        if isAdaptiveStateCovariance && ~usePerfectFrozenTruthState
            adaptiveStateCovarianceError = ...
                innovation - jacobian * (x_k_k - x_k_k_1);
            adaptiveStateCovarianceErrorMemory = ...
                [adaptiveStateCovarianceError ...
                adaptiveStateCovarianceErrorMemory(:, 1:end-1)];
            if k > adaptiveStateCovarianceWindow
                innovationErrorCovariance = ...
                    (adaptiveStateCovarianceErrorMemory * ...
                    adaptiveStateCovarianceErrorMemory') / ...
                    adaptiveStateCovarianceWindow;
                Q = K_k * innovationErrorCovariance * K_k';
            end
        end

        P_k_k = (eye(stateDimension) - K_k * jacobian) * P_k_k_1;
        if ~useAcausalTaps
            x_k_k(acausalTapStateIndices) = 0;
            P_k_k(acausalTapStateIndices, :) = 0;
            P_k_k(:, acausalTapStateIndices) = 0;
        end
    else
        x_k_k = x_k_k_1;
        P_k_k = P_k_k_1;
    end

    %% EKF Projection Ahead Step
    x_k_k_1 = x_k_k;
    x_k_k_1(delayStateIndex) = real(x_k_k_1(delayStateIndex));
    if ~useAcausalTaps
        x_k_k_1(acausalTapStateIndices) = 0;
    end
    P_k_k_1 = P_k_k + Q;
    if ~useAcausalTaps
        P_k_k_1(acausalTapStateIndices, :) = 0;
        P_k_k_1(:, acausalTapStateIndices) = 0;
    end

    delayEstimateRecord(:, k) = x_k_k(delayStateIndex);
    channelStateRecord(:, k) = x_k_k(tapStateIndices);
end

%% Plots
lineWidth = 2;
fontSize = 13;

epochVector = 1:simulationSteps;
timeMs = (epochVector - 1) * epoch * 1e3;
chipPeriod = 1 / configuration.chippingFrequency;
truthSampleIndex = round(epochVector * samplesTotal);
trueDelayRecord = LOSDelay(truthSampleIndex).';
trueDelayEstimationErrorRecord = trueDelayRecord - delayEstimateRecord;

forwardTapIndices = (q + 1):tapCount;
numberOfForwardTaps = numel(forwardTapIndices);
forwardTapNumbers = 0:(numberOfForwardTaps - 1);
forwardTapDelayTc = forwardTapNumbers / samplesPerChip;
estimatedForwardChannel = channelStateRecord(forwardTapIndices, :);
estimatedForwardChannelMean = mean(estimatedForwardChannel, 2);
estimatedForwardChannelRealStd = std(real(estimatedForwardChannel), 0, 2);
estimatedForwardChannelImagStd = std(imag(estimatedForwardChannel), 0, 2);
trueForwardChannel = truthChannelState(forwardTapIndices);
channelPlotRows = ceil(sqrt(numberOfForwardTaps));
channelPlotColumns = ceil(numberOfForwardTaps / channelPlotRows);

innovationStdRecord = std(innovationRecord(1:C, :), 1, 1);
figure(Name="STD of the innovations", NumberTitle="off");
hold on;
plot(timeMs, innovationStdRecord, 'LineWidth', lineWidth);
plot(timeMs, zeros(1, simulationSteps), '--', 'LineWidth', lineWidth);
legend({"Innovation STD", "Truth"});
ylabel("Standard deviation of the innovations");
xlabel("Time [ms]");
hold off;

figure(Name="Middle tap of the innovation sequence", NumberTitle="off");
hold on;
plot(timeMs, real(innovationRecord(middleSample, :)), 'LineWidth', lineWidth);
plot(timeMs, imag(innovationRecord(middleSample, :)), 'LineWidth', lineWidth);
plot(timeMs, zeros(1, simulationSteps), '--', 'LineWidth', lineWidth);
legend({"Real", "Imaginary", "Truth"});
ylabel("Innovation sequence of the middle tap");
xlabel("Time [ms]");
hold off;

figure(Name="Delay Estimation", NumberTitle="off");
hold on;
plot(timeMs, (delayEstimateRecord - trueDelayRecord(1)) / chipPeriod, ...
    'LineWidth', lineWidth);
plot(timeMs, (trueDelayRecord - trueDelayRecord(1)) / chipPeriod, ...
    '--', 'LineWidth', lineWidth);
legend({"EKF's estimated delay", "True delay"});
ylabel("LOS Delay [T_c]");
xlabel("Time [ms]");
hold off;
%
% figure(Name="Real Kalman Gain Elements", NumberTitle="off");
% hold on;
% for i = 1:size(kalmanGainRecord, 1)
%     for j = 1:size(kalmanGainRecord, 2)
%         plot(timeMs, squeeze(real(kalmanGainRecord(i, j, :))), ...
%             'DisplayName', sprintf('K_{%d,%d}', i, j), ...
%             'LineWidth', lineWidth);
%     end
% end
% legend show;
% ylabel("Real Kalman Gain Elements");
% xlabel("Time [ms]");
% set(gca, "FontSize", fontSize);
% hold off;
%
% figure(Name="Imaginary Kalman Gain Elements", NumberTitle="off");
% hold on;
% for i = 1:size(kalmanGainRecord, 1)
%     for j = 1:size(kalmanGainRecord, 2)
%         plot(timeMs, squeeze(imag(kalmanGainRecord(i, j, :))), ...
%             'DisplayName', sprintf('K_{%d,%d}', i, j), ...
%             'LineWidth', lineWidth);
%     end
% end
% legend show;
% ylabel("Imaginary Kalman Gain Elements");
% xlabel("Time [ms]");
% set(gca, "FontSize", fontSize);
% hold off;

% figure(Name="secondary taps", NumberTitle="off");
% hold on;
% for i = 1:q
%     tapIndex = i + q + 1;
%     plot(timeMs, abs(channelStateRecord(tapIndex, :)), ...
%         'LineWidth', lineWidth, ...
%         'DisplayName', sprintf("Estimate tap %+d", i));
%     plot(timeMs, abs(truthChannelState(tapIndex)) * ...
%         ones(1, simulationSteps), '--', ...
%         'LineWidth', lineWidth, ...
%         'DisplayName', sprintf("Truth tap %+d", i));
% end
% legend show;
% ylabel("Magnitude");
% xlabel("Time [ms]");
% hold off;
%
% figure(Name="main tap", NumberTitle="off");
% hold on;
% plot(timeMs, real(channelStateRecord(q+1, :)), 'LineWidth', lineWidth);
% plot(timeMs, real(truthChannelState(q+1)) * ...
%     ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
% legend({"Estimate", "Truth"});
% ylabel("Real");
% xlabel("Time [ms]");
% hold off;

figure(Name="Real Channel Tap Estimates", NumberTitle="off");
tiledlayout(channelPlotRows, channelPlotColumns);
for tapIndex = forwardTapIndices
    nexttile;
    hold on;
    plot(timeMs, real(channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(timeMs, real(truthChannelState(tapIndex)) * ...
        ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
    hold off;
    title(sprintf("Tap %+d", tapIndex - q - 1));
    ylabel("Real");
    xlabel("Time [ms]");
    set(gca, "FontSize", fontSize);
end
legend({"Estimate", "Truth"});

figure(Name="Imaginary Channel Tap Estimates", NumberTitle="off");
tiledlayout(channelPlotRows, channelPlotColumns);
for tapIndex = forwardTapIndices
    nexttile;
    hold on;
    plot(timeMs, imag(channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(timeMs, imag(truthChannelState(tapIndex)) * ...
        ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
    hold off;
    title(sprintf("Tap %+d", tapIndex - q - 1));
    ylabel("Imag");
    xlabel("Time [ms]");
    set(gca, "FontSize", fontSize);
end
legend({"Estimate", "Truth"});

figure(Name="Channel Coefficient Summary", NumberTitle="off");
tiledlayout(1, 2);

nexttile;
realBar = bar(forwardTapNumbers, ...
    [real(estimatedForwardChannelMean) real(trueForwardChannel)]);
realBar(1).FaceColor = [0.55 0.75 0.95];
realBar(2).FaceColor = [0 0.5 0];
hold on;
errorbar(realBar(1).XEndPoints, real(estimatedForwardChannelMean), ...
    estimatedForwardChannelRealStd, 'LineStyle', 'none', ...
    'Color', [0 0.4470 0.7410], 'LineWidth', lineWidth);
hold off;
grid on;
legend({"Estimated Channel Coefficients", ...
    "True Channel Coefficients"});
title("Real Channel Coefficients");
ylabel("Channel Coefficients");
xticks(forwardTapNumbers);
xticklabels(compose("h_{%d}", forwardTapNumbers));
set(gca, "FontSize", fontSize);

nexttile;
imagBar = bar(forwardTapNumbers, ...
    [imag(estimatedForwardChannelMean) imag(trueForwardChannel)]);
imagBar(1).FaceColor = [0.55 0.75 0.95];
imagBar(2).FaceColor = [0 0.5 0];
hold on;
errorbar(imagBar(1).XEndPoints, imag(estimatedForwardChannelMean), ...
    estimatedForwardChannelImagStd, 'LineStyle', 'none', ...
    'Color', [0 0.4470 0.7410], 'LineWidth', lineWidth);
hold off;
grid on;
legend({"Estimated Channel Coefficients", ...
    "True Channel Coefficients"});
title("Imaginary Channel Coefficients");
ylabel("Channel Coefficients");
xticks(forwardTapNumbers);
xticklabels(compose("h_{%d}", forwardTapNumbers));
set(gca, "FontSize", fontSize);

figure(Name="Channel Impulse Response History", NumberTitle="off");
[timeGrid, tapDelayGrid] = meshgrid(timeMs, forwardTapDelayTc);
tiledlayout(1, 2);

nexttile;
hold on;
surf(timeGrid, tapDelayGrid, real(estimatedForwardChannel), ...
    'EdgeColor', 'none');
surf(timeGrid, tapDelayGrid, ...
    repmat(real(trueForwardChannel), 1, simulationSteps), ...
    'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.85, ...
    'EdgeColor', 'none');
hold off;
grid on;
colormap(turbo);
colorbar;
view(42, 28);
camproj("perspective");
axis tight;
legend({"Estimate", "Truth"});
title("Real Channel Impulse Response");
xlabel("Time [ms]");
ylabel("Delay [T_c]");
zlabel("Real");
set(gca, "FontSize", fontSize);

nexttile;
hold on;
surf(timeGrid, tapDelayGrid, imag(estimatedForwardChannel), ...
    'EdgeColor', 'none');
surf(timeGrid, tapDelayGrid, ...
    repmat(imag(trueForwardChannel), 1, simulationSteps), ...
    'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.85, ...
    'EdgeColor', 'none');
hold off;
grid on;
colormap(turbo);
colorbar;
view(42, 28);
camproj("perspective");
axis tight;
legend({"Estimate", "Truth"});
title("Imaginary Channel Impulse Response");
xlabel("Time [ms]");
ylabel("Delay [T_c]");
zlabel("Imaginary");
set(gca, "FontSize", fontSize);

figure(Name="Final Channel Tap Estimates", NumberTitle="off");
hold on;
plot(real(truthChannelState(forwardTapIndices)), ...
    imag(truthChannelState(forwardTapIndices)), 'x', ...
    'LineWidth', lineWidth, 'MarkerSize', 10);
plot(real(channelStateRecord(forwardTapIndices, end)), ...
    imag(channelStateRecord(forwardTapIndices, end)), ...
    'o', 'LineWidth', lineWidth, 'MarkerSize', 8);
for tapIndex = forwardTapIndices
    text(mean(real(channelStateRecord(tapIndex))), ...
        mean(imag(channelStateRecord(tapIndex))), ...
        sprintf(" %+d", tapIndex - q - 1));
end
grid on;
axis equal;
legend({"Truth", "Mean estimate"});
xlabel("Real");
ylabel("Imaginary");
set(gca, "FontSize", fontSize);
hold off;

% VariationRecord = zeros(size(kalmanGainRecord, 1), simulationSteps);
% for i = 1:simulationSteps
%     VariationRecord(:, i) = kalmanGainRecord(:, :, i) * innovationRecord(:, i);
% end

function shiftedCorrelations = getShiftedCorrelationsFromBank( ...
    correlatorBank, delay, configuration, channelTapDelays)

samplesPerEpoch = size(correlatorBank, 2);
channelReplicas = zeros(numel(channelTapDelays), samplesPerEpoch);
for col = 1:numel(channelTapDelays)
    channelReplicas(col, :) = getCodeReplica( ...
        configuration, delay + channelTapDelays(col)).';
end

shiftedCorrelations = correlatorBank * channelReplicas.';

end
