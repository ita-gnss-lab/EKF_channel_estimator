clearvars; clc; close all;

addpath(genpath(fullfile("..", "..","EKF_channel_estimator")));

load config_no_doppler.mat
rng(26437226); 
samplesPerChip = 16;
configuration.samplingFrequency = samplesPerChip * configuration.chippingFrequency;

%% Stress-test truth
configuration.carrierToNoiseDensityRatio = 42;
configuration.dopplerProfile(1:3) = [0 0 0];
configuration.dopplerProfile(1) = ...
    -2*pi*configuration.carrierFrequency*1e-4;
numberOfCausalTruthTaps = 19;
truthTapOrder = 0:(numberOfCausalTruthTaps - 1);
configuration.tdl_channel = 0.8 * ...
    exp((-0.85 + 1j*pi/5) * truthTapOrder);
%% Parameters
simulationSteps = 1500;
useTapEnergyConstraint = false;
useAcausalTaps = true;
usePerfectFrozenTruthState = false;
constraint_noise = 10^(-2.63);%10.^[-2.63 -3.44 -4 -4.28 -4.49 -4.63 -4.85 -5.2 -5.37 -5.88];
q = numberOfCausalTruthTaps;
C = 2*q + 1;
middleSample = q + 1;
acausalTapStateIndices = 5:(4+q);
acausalTapChannelIndices = 1:q;
epoch = configuration.totalChips / configuration.chippingFrequency;
configuration.correlatorHalfSpan = q;
samplesTotal = epoch*configuration.samplingFrequency;
timeSupport = (0:(samplesTotal - 1)).' * (1/configuration.samplingFrequency);
beta = -1 / (2 * pi * configuration.carrierFrequency);
WienerStatesSelection = 1:4;
trueDelay = -configuration.dopplerProfile(1) / ...
    (2*pi*configuration.carrierFrequency);
initialDelayErrorSamples = 0;
initialDelayEstimate = trueDelay + ...
    initialDelayErrorSamples / configuration.samplingFrequency;
initialDelayError = trueDelay - initialDelayEstimate;
delayErrorStdSamples = max(2, abs(initialDelayErrorSamples));

truthChannelAcausal = zeros(2*q + 1, 1);
numberOfTruthTaps = min(numel(configuration.tdl_channel), q + 1);
truthChannelAcausal(q + 1:q + numberOfTruthTaps) = ...
    configuration.tdl_channel(1:numberOfTruthTaps).';
truthChannelState = truthChannelAcausal;

%% Covariances 
% Convert CN0 from dB-Hz to linear scale
carrierToNoiseRatioLinear = 10^(configuration.carrierToNoiseDensityRatio / 10);
% Compute the noise variance
thermalNoiseVarianceSquared = configuration.samplingFrequency / carrierToNoiseRatioLinear;
sigma2Vec = [1e-1 1e1 0 0 1e-4];
Q = getStateCovarianceMatrix_acausal(...
    sigma2Vec, ...
    epoch, ...
    beta, ...
    q);
if ~useAcausalTaps
    Q(acausalTapStateIndices, :) = 0;
    Q(:, acausalTapStateIndices) = 0;
end
if usePerfectFrozenTruthState
    Q(:) = 0;
end
correlatorBank = buildCorrelatorBank(configuration, 0, q);

R  = (thermalNoiseVarianceSquared / samplesTotal.^2) * ...
        (correlatorBank * correlatorBank.');
if useTapEnergyConstraint
    R = [R zeros(2*q+1, 1); zeros(1, 2*q+1) constraint_noise];
end

%% State History Vectors
LQGStateRecord = zeros(4, simulationSteps);
errorStateRecord = zeros(4, simulationSteps);
channelStateRecord = zeros(2*q + 1, simulationSteps);
measurementDimension = C + double(useTapEnergyConstraint);
innovationRecord = zeros(measurementDimension, simulationSteps);
kalmanGainRecord = zeros(4 + 2*q + 1, measurementDimension, simulationSteps);
constraintRecord = zeros(1, simulationSteps);

%% Cost Functions
% I4 = eye(4);
% 
% % Bryson-style targets (tune as needed)
% sig_tau  = 1e-4;    % s
% sig_phi  = 2e-3;    % rad
% sig_nu   = 1e-2;    % rad/
% sig_nud  = 2e-1;    % rad/s^2
% 
% T_e_0 = diag([1/sig_tau^2, 1/sig_phi^2, 1/sig_nu^2, 1/sig_nud^2]);
% T_u_0 = 1e-2*eye(4);  % small => fast
% 
% kappa = 5; % speed knob
% T_e = kappa*T_e_0;
% T_u = T_u_0/kappa;

relation = 10;
delayCostScale = abs(beta);
T_e = relation * diag([delayCostScale, 1, 1/epoch, 1/epoch^2]);
T_u = diag([delayCostScale, 1, 1/epoch, 1/epoch^2]);

%% Transition Matrices
[F_W, F_H] = getModelTransitionMatrix_acausal(epoch, q, beta);

%% Coupling Matrix for Control Signal
B_LQG = eye(4);

%% IDARE Solution
[~, L, ~] = idare(F_W, B_LQG, T_e, T_u, [], []);

%% Full transition matrix
F = blkdiag(F_W, F_H);


%% Initialization
% NOTE: I changed from x_k_k to x_k_k_1, because, in fact the
% initialization uses x[1|0].
x_k_k_1 = zeros(4 + (2*q+1), 1);
main_tap = false(size(x_k_k_1));
main_tap(4 + q + 1) = true;
x_k_k_1(main_tap) = 1;
other_taps = true(size(x_k_k_1));
other_taps(1:4) = false;
other_taps(main_tap) = false;
x_k_k_1(other_taps) = 0;
if usePerfectFrozenTruthState
    x_k_k_1(1:4) = 0;
    x_k_k_1(5:end) = truthChannelState;
end

% NOTE: I changed from x_k_k to P_k_k_1, because, in fact the
% initialization uses P[1|0].  1e-1, 0, (50)^2/12, (0.1)^2/12,
delayErrorStd = delayErrorStdSamples / configuration.samplingFrequency;
phaseErrorStd = 0;
channelCovarianceMatrix = 1e-2 * eye(2*q + 1);
channelCovarianceMatrix(q + 1, q + 1) = 1e-1;
if usePerfectFrozenTruthState
    channelCovarianceMatrix(:) = 0;
end
if ~useAcausalTaps
    channelCovarianceMatrix(acausalTapChannelIndices, :) = 0;
    channelCovarianceMatrix(:, acausalTapChannelIndices) = 0;
end
P_k_k_1 = blkdiag( ...
    delayErrorStd^2, ...
    10*phaseErrorStd^2, ...
    0, ...
    0, ...
    channelCovarianceMatrix);
% P_k_k_1 = blkdiag(1e-1, 0, 0, 0, zeros(1 + q));

initialPhaseError = phaseErrorStd;
x_LQG_k = [initialDelayEstimate, ...
    configuration.dopplerProfile(1) + initialPhaseError, ...
    2*pi*configuration.dopplerProfile(2), ...
    2*pi*configuration.dopplerProfile(3)].';

u_LQG = L * x_k_k_1(WienerStatesSelection);
 u_LQG(2) = 0;
u_LQG(3:4) = 0;

%% Simulate Signal
configuration.addNoise = false;
[simulatedSignal, ~, LOSPhase, LOSDelay] = gnssReceivedSignal(configuration, simulationSteps + 1);

%% Simulation
plotMeasures = false;
correlatorTaps = -q:1:q;
% NOTE: (Rodrigo): Changed the main loop to match algorithm 1 of my report.
for k = 1 : simulationSteps
    %% Signal 
    receivedSignal = simulatedSignal(((k - 1) * samplesTotal + 1: k * samplesTotal));

    %% LQG Controller
    % Update LQG state 
    x_LQG_k = F_W * x_LQG_k + B_LQG * u_LQG;
    % x_LQG_k(2) = configuration.dopplerProfile(1);
    % x_LQG_k(3:4) = [
    %     2*pi*configuration.dopplerProfile(2);
    %     2*pi*configuration.dopplerProfile(3)];

    LQGStateRecord(:, k) = x_LQG_k(1:4);

    % Carrier Wipe-Off
    phi_T = x_LQG_k(2) + x_LQG_k(3) * timeSupport + 0.5 * x_LQG_k(4) * timeSupport.^2;
    d_k = exp(1j * phi_T);
    wipedSignal = receivedSignal .* conj(d_k);

    %% Kalman filter
    if k > 1
        % EKF's Update Step
        correlatorBank = buildCorrelatorBank(configuration, x_LQG_k(1), q);
        z_k = correlatorBank * wipedSignal / samplesTotal;
        constraint_value = sum(abs(x_k_k_1(other_taps)).^2)/((q-1)*abs(x_k_k_1(main_tap))^2);
        z_hat_k_aux = measurementFunction_acausal( ...
            x_k_k_1, configuration, q, x_LQG_k(1)) / samplesTotal;
        z_hat_k = z_hat_k_aux;
        if useTapEnergyConstraint
            z_k = [z_k; 0];
            z_hat_k = [z_hat_k; constraint_value];
        end
        constraintRecord(:, k) = constraint_value;
        
        if plotMeasures
            % ---- Plot routine -----
            plot(correlatorTaps, real(z_k(1:C)));
            hold on;
            plot(correlatorTaps, real(z_hat_k(1:C)));
            plot(correlatorTaps, imag(z_k(1:C)));
            plot(correlatorTaps, imag(z_hat_k(1:C)));
            plot(correlatorTaps(q+1:end), real(z_k(q+1:C)), 'o');
            plot(correlatorTaps(q+1:end), real(z_hat_k(q+1:C)), 'x');
            hold off;
            ylabel('Real and Imag parts of z_k and z_k_hat');
            xlabel('Correlator tap');
            legend({'Real $z[k]$', 'Real $\hat{z}[k]$', ...
                'Imag $z[k]$', 'Imag $\hat{z}[k]$', ...
                'Causal real $z[k]$', 'Causal real $\hat{z}[k]$'}, ...
                'Interpreter','latex');
            pause(0.01);
        end
        
        % Compute Jacobian
        % Delay term now follows Φ_pp(ετ + (l - m)Ts) as in the analytical model.
        delayJacobian = delayJacobianFunctionSimplified_acausal( ...
            x_k_k_1(1), ...
            x_k_k_1(5:end), ...
            1 / configuration.samplingFrequency, ...
            q, ...
            1 / configuration.chippingFrequency, ...
            exp(1j * x_k_k_1(2)) ...
        );
        phaseJacobian = 1j * z_hat_k_aux;
        dopplerJacobian = zeros(2*q + 1, 2);
        channelOrder = (numel(x_k_k_1(5:end)) - 1)/2;
        channelWeightsJacobian = exp(1j * x_k_k_1(2)) .* ...
            getShiftedCorrelations_acausal( ...
                x_k_k_1(1), q, configuration, channelOrder, x_LQG_k(1)) / ...
            samplesTotal;
        % todo - add the jacobian of the constraint
        LOSParcel = -x_k_k_1(main_tap)*sum(abs(x_k_k_1(other_taps)).^2)/(abs(x_k_k_1(main_tap))^2);
        tapsParcel = (2/((q-1)*abs(x_k_k_1(main_tap))^2))*[x_k_k_1(4+1:4+q);LOSParcel;x_k_k_1(4+q+2:end)];
        constraintLine = [0 0 0 0 tapsParcel'];
        jacobian = [delayJacobian ...
            phaseJacobian ...
            dopplerJacobian ...
            channelWeightsJacobian];
        if useTapEnergyConstraint
            jacobian = [jacobian; constraintLine];
        end
        
        % Compute Kalman Gain
        K_k = P_k_k_1 * jacobian'...
            *((jacobian * P_k_k_1 * jacobian' + R) \ eye(measurementDimension));
        kalmanGainRecord(:,:, k) = K_k;

        % Obtain the innovation
        innovation = z_k - z_hat_k;
        innovationRecord(:, k) = innovation;

        % EKF's state update
        x_k_k = x_k_k_1 + K_k * innovation;
        x_k_k(WienerStatesSelection) = real(x_k_k(WienerStatesSelection));
        % x_k_k(2) = 0;
        x_k_k(3:4) = 0;
        if ~useAcausalTaps
            x_k_k(acausalTapStateIndices) = 0;
        end
 
        % EKF's covariance matrix update
        P_k_k = (eye(2*q + 1 + 4) - K_k*jacobian) * P_k_k_1;
        if ~useAcausalTaps
            P_k_k(acausalTapStateIndices, :) = 0;
            P_k_k(:, acausalTapStateIndices) = 0;
        end
        
        % LQG control vector computation
        u_LQG = L * x_k_k(WienerStatesSelection);
        % u_LQG(2) = 0;
        u_LQG(3:4) = 0;
    else
        % Initialization procedure
        x_k_k = x_k_k_1;
        P_k_k = P_k_k_1;
    end

    % EKF's Projection Ahead Step
    x_k_k_1 = F * x_k_k;  
    x_k_k_1(WienerStatesSelection) = real(x_k_k_1(WienerStatesSelection));
    % x_k_k_1(2) = 0;
    x_k_k_1(3:4) = 0;
    if ~useAcausalTaps
        x_k_k_1(acausalTapStateIndices) = 0;
    end
    P_k_k_1 = F * P_k_k * F' + Q;
    if ~useAcausalTaps
        P_k_k_1(acausalTapStateIndices, :) = 0;
        P_k_k_1(:, acausalTapStateIndices) = 0;
    end

    errorStateRecord(:, k) = x_k_k(1:4);
    channelStateRecord(:, k) = x_k_k(5:end);

end
%% Plots

lineWidth = 2;
fontSize = 13;

epochVector = 1:simulationSteps;
truthSampleIndex = round(epochVector * samplesTotal);
trueDelayRecord = LOSDelay(truthSampleIndex).';
truePhaseRecord = LOSPhase(truthSampleIndex).';
trueDopplerRecord = 2*pi*configuration.dopplerProfile(2) * ...
    ones(1, simulationSteps);
trueDelayErrorRecord = trueDelayRecord - LQGStateRecord(1, :);
truePhaseErrorRecord = truePhaseRecord - LQGStateRecord(2, :);
trueDopplerErrorRecord = trueDopplerRecord - LQGStateRecord(3, :);

truthChannelAcausal = zeros(2*q + 1, 1);
numberOfTruthTaps = min(numel(configuration.tdl_channel), q + 1);
truthChannelAcausal(q + 1:q + numberOfTruthTaps) = ...
    configuration.tdl_channel(1:numberOfTruthTaps).';
truthChannelState = truthChannelAcausal;
forwardTapIndices = (q + 1):(2*q + 1);
numberOfForwardTaps = numel(forwardTapIndices);
channelPlotRows = ceil(sqrt(numberOfForwardTaps));
channelPlotColumns = ceil(numberOfForwardTaps / channelPlotRows);

%Observe the STD of the innovation sequence time series\
innovationStdRecord = std(innovationRecord(1:C,:),1,1);
figure(Name="STD of the innovations", NumberTitle="off");
hold on;
plot(epochVector, innovationStdRecord, 'LineWidth', lineWidth);
plot(epochVector, zeros(1, simulationSteps), '--', 'LineWidth', lineWidth);
legend({"Innovation STD", "Truth"});
ylabel("Standard deviation of the innovations");
xlabel("Epochs (Simulation Steps)");
hold off;

% Observe the innovation sequence time series
figure(Name="Middle tap of the innovation sequence", NumberTitle="off");
hold on;
plot(epochVector, real(innovationRecord(middleSample,:)), 'LineWidth', lineWidth);
plot(epochVector, imag(innovationRecord(middleSample,:)), 'LineWidth', lineWidth);
plot(epochVector, zeros(1, simulationSteps), '--', 'LineWidth', lineWidth);
legend({"Real", "Imaginary", "Truth"});
ylabel("Innovation sequence of the middle tap");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="Delay Estimation", NumberTitle="off");
hold on;
plot(epochVector, LQGStateRecord(1,:), 'LineWidth', lineWidth);
plot(epochVector, trueDelayRecord, '--', 'LineWidth', lineWidth);
legend({"LQG's estimated delay", "True delay"});
ylabel("Delay estimate");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="Delay Error State", NumberTitle="off");
hold on;
plot(errorStateRecord(1, :), 'LineWidth', lineWidth);
plot(epochVector, trueDelayErrorRecord, '--', 'LineWidth', lineWidth);
plot(epochVector, zeros(1, simulationSteps), ':', 'LineWidth', lineWidth);
legend({"EKF's estimated delay error", "True delay error", "Zero line"});
ylabel("Delay error estimate");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="Phase Estimation", NumberTitle="off");
hold on;
plot(epochVector, LQGStateRecord(2,:), 'LineWidth', lineWidth);
plot(epochVector, truePhaseRecord, '--', 'LineWidth', lineWidth);
legend({"LQG's estimated phase", "True Phase"});
ylabel("Phase estimate");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="Phase Error State", NumberTitle="off");
hold on;
plot(errorStateRecord(2, :), 'LineWidth', lineWidth);
plot(epochVector, truePhaseErrorRecord, '--', 'LineWidth', lineWidth);
plot(epochVector, zeros(1, simulationSteps), ':', 'LineWidth', lineWidth);
legend({"EKF's estimated phase error", "True phase error", "Zero line"});
ylabel("Phase error estimate");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="Doppler Estimation", NumberTitle="off");
hold on;
plot(epochVector, LQGStateRecord(3,:), 'LineWidth', lineWidth);
plot(epochVector, trueDopplerRecord, '--', 'LineWidth', lineWidth);
legend({"LQG's estimated Doppler frequency", "True Doppler frequency"});
ylabel("Doppler estimate");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="Doppler Error State", NumberTitle="off");
hold on;
plot(errorStateRecord(3, :), 'LineWidth', lineWidth);
plot(epochVector, trueDopplerErrorRecord, '--', 'LineWidth', lineWidth);
plot(epochVector, zeros(1, simulationSteps), ':', 'LineWidth', lineWidth);
ylabel("Doppler error estimate");
xlabel("Epochs (Simulation Steps)");
legend({"EKF's estimated Doppler error", "True Doppler error", "Zero line"});
hold off;

figure(Name="Real Kalman Gain Elements", NumberTitle="off");
hold on;
for i = 1:size(kalmanGainRecord, 1)
    for j = 1:size(kalmanGainRecord, 2)
        plot(epochVector, squeeze(real(kalmanGainRecord(i, j, :))), 'DisplayName', sprintf('K_{%d,%d}', i, j), 'LineWidth', lineWidth);
    end
end
legend show;
ylabel("Real Kalman Gain Elements");
xlabel("Epochs (Simulation Steps)");
set(gca, "FontSize", fontSize);
hold off;

figure(Name="Imaginary Kalman Gain Elements", NumberTitle="off");
hold on;
for i = 1:size(kalmanGainRecord, 1)
    for j = 1:size(kalmanGainRecord, 2)
        plot(epochVector, squeeze(imag(kalmanGainRecord(i, j, :))), 'DisplayName', sprintf('K_{%d,%d}', i, j), 'LineWidth', lineWidth);
    end
end
legend show;
ylabel("Imaginary Kalman Gain Elements");
xlabel("Epochs (Simulation Steps)");
set(gca, "FontSize", fontSize);
hold off;

% figure(Name="Channel Weights Over Time", NumberTitle="off");
% hold on;
% for i = 1:3
%     plot(abs(channelStateRecord(i, :)), 'LineWidth', lineWidth);
% end
% yyaxis right
% plot(constraintRecord, 'LineWidth', lineWidth);
% ylabel("Constraint");
% xlabel("Epochs (Simulation Steps)");
% hold off;

figure(Name="secondary taps", NumberTitle="off");
hold on;
for i = 1:q
    tapIndex = i + q + 1;
    plot(abs(channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth, ...
        'DisplayName', sprintf("Estimate tap %+d", i));
    plot(epochVector, abs(truthChannelState(tapIndex)) * ...
        ones(1, simulationSteps), '--', ...
        'LineWidth', lineWidth, ...
        'DisplayName', sprintf("Truth tap %+d", i));
end
legend show;
ylabel("Magnitude");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="main tap", NumberTitle="off");
hold on; 
plot(real(channelStateRecord(q+1, :)), 'LineWidth', lineWidth);
plot(epochVector, real(truthChannelState(q+1)) * ...
    ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
legend({"Estimate", "Truth"});
ylabel("Real");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="Real Channel Tap Estimates", NumberTitle="off");
tiledlayout(channelPlotRows, channelPlotColumns);
for tapIndex = forwardTapIndices
    nexttile;
    hold on;
    plot(epochVector, real(channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, real(truthChannelState(tapIndex)) * ...
        ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
    hold off;
    title(sprintf("Tap %+d", tapIndex - q - 1));
    ylabel("Real");
    xlabel("Epochs");
    set(gca, "FontSize", fontSize);
end
legend({"Estimate", "Truth"});

figure(Name="Imaginary Channel Tap Estimates", NumberTitle="off");
tiledlayout(channelPlotRows, channelPlotColumns);
for tapIndex = forwardTapIndices
    nexttile;
    hold on;
    plot(epochVector, imag(channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, imag(truthChannelState(tapIndex)) * ...
        ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
    hold off;
    title(sprintf("Tap %+d", tapIndex - q - 1));
    ylabel("Imag");
    xlabel("Epochs");
    set(gca, "FontSize", fontSize);
end
legend({"Estimate", "Truth"});

figure(Name="Final Channel Tap Estimates", NumberTitle="off");
hold on;
plot(real(truthChannelState(forwardTapIndices)), ...
    imag(truthChannelState(forwardTapIndices)), 'x', ...
    'LineWidth', lineWidth, 'MarkerSize', 10);
plot(real(channelStateRecord(forwardTapIndices, end)), ...
    imag(channelStateRecord(forwardTapIndices, end)), ...
    'o', 'LineWidth', lineWidth, 'MarkerSize', 8);
for tapIndex = forwardTapIndices
    text(real(channelStateRecord(tapIndex, end)), ...
        imag(channelStateRecord(tapIndex, end)), ...
        sprintf(" %+d", tapIndex - q - 1));
end
grid on;
axis equal;
legend({"Truth", "Final estimate"});
xlabel("Real");
ylabel("Imaginary");
set(gca, "FontSize", fontSize);
hold off;

VariationRecord = zeros(size(kalmanGainRecord, 1), simulationSteps);
for i = 1:simulationSteps
    VariationRecord(:,i) = kalmanGainRecord(:,:,i)*innovationRecord(:,i);
end
