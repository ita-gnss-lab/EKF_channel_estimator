clearvars; clc; close all;

addpath(genpath(fullfile("..", "..","EKF_channel_estimator")));

load config_cte_doppler.mat
rng(26437226);

%% Parameters
simulationSteps = 500;
q = 4;
C = 2*q + 1;
middleSample = q + 1;
epoch = configuration.totalChips / configuration.chippingFrequency;
configuration.correlatorHalfSpan = q;
samplesTotal = epoch*configuration.samplingFrequency;
timeSupport = (0:(samplesTotal - 1)).' * (1/configuration.samplingFrequency);
beta = -1 / (2 * pi * configuration.carrierFrequency);
WienerStatesSelection = 1:4;

%% Covariances 
% Convert CN0 from dB-Hz to linear scale
carrierToNoiseRatioLinear = 10^(configuration.carrierToNoiseDensityRatio / 10);
% Compute the noise variance
thermalNoiseVarianceSquared = configuration.samplingFrequency / carrierToNoiseRatioLinear;

sigma2Vec = [1e-1 1e-1 1e-2 1e-3 1e-6];
Q = getStateCovarianceMatrix(...
    sigma2Vec, ...
    epoch, ...
    beta, ...
    q);
correlatorBank = buildCorrelatorBank(configuration, 0, q);
R = (thermalNoiseVarianceSquared / samplesTotal.^2) * ...
        (correlatorBank * correlatorBank.');

%% State History Vectors
LQGStateRecord = zeros(4, simulationSteps);
errorStateRecord = zeros(4, simulationSteps);
innovationRecord = zeros(C, simulationSteps);
kalmanGainRecord = zeros(4 + q + 1,C, simulationSteps);
channelStateRecord = zeros(q + 1, simulationSteps);

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
T_e = relation * diag([beta, 1, 1/epoch, 1/epoch^2]);
T_u = diag([beta, 1, 1/epoch, 1/epoch^2]);

%% Transition Matrices
[F_W, F_H] = getModelTransitionMatrix(epoch, q, beta);

%% Coupling Matrix for Control Signal
B_LQG = eye(4);

%% IDARE Solution
[~, L, ~] = idare(F_W, B_LQG, T_e, T_u, [], []);

%% Full transition matrix
F = blkdiag(F_W, F_H);

%% --- Random channel taps: Rician with exponential PDP ---
rng(42);                               % reproducible
K_dB = 25;
K = 10^(K_dB/10);
tau = 1.2; % PDP decay (samples)
totalPow = 1; % total tap power (sum |h|^2)

pdp = exp(-(0:q + 1-1)/tau);                % exponential PDP
pdp = pdp / sum(pdp) * totalPow;        % normalize to target total power

losIdx = 1; % LOS on tap 0 (use 1-based idx); choose another if you want
h_los = zeros(q + 1,1);
h_los(losIdx) = sqrt(K/(K+1) * pdp(losIdx));

% scattered (Rayleigh) on all taps
h_scatt = (randn(q + 1,1) + 1j*randn(q + 1,1))/sqrt(2) .* sqrt((1/(K+1)) * pdp(:));

channelTaps = (h_los + h_scatt).';      % row vector to match your code

%% Initialization
% NOTE: I changed from x_k_k to x_k_k_1, because, in fact the
% initialization uses x[1|0].
x_k_k_1 = zeros(4 + q + 1, 1);

% EKF channel state init: noisy around truth
sigmaInit = 1e-1; % relative RMS vs PDP (tune)

% scale noise per PDP so early taps tend to have larger prior variance
pdp_est = abs(channelTaps).^2;
pdp_est = pdp_est / max(pdp_est);  % simple shape; avoid huge values

x_k_k_1(5:5 + q + 1 - 1) = channelTaps(:) + ...
    sigmaInit * ((randn(q + 1,1)+1j*randn(q + 1,1))/sqrt(2)) .* sqrt(pdp_est(:));

% HACK: I zeroed this initial covariance matrix to my analysis about the
% phase estimation.
channelCovarianceMatrix = 1e-8 * eye(q + 1); %0.000001 * eye(1 + q);

% NOTE: I changed from x_k_k to P_k_k_1, because, in fact the
% initialization uses P[1|0].
P_k_k_1 = blkdiag(1e-1, 0, (50)^2/12, (0.1)^2/12, channelCovarianceMatrix); 
% P_k_k_1 = blkdiag(1e-1, 0, 0, 0, zeros(1 + q));

phaseError = 0.5;
DopplerError = 25;
x_LQG_k = [1.005e-4, ...
    configuration.dopplerProfile(1) + phaseError, ...
    2*pi*(configuration.dopplerProfile(2) + DopplerError), ...
    2*pi*configuration.dopplerProfile(3)].';

u_LQG = L * x_k_k_1(WienerStatesSelection);

%% Simulate Signal
configuration.addNoise = true;
[simulatedSignalRaw, ~, LOSPhase, LOSDelay] = gnssReceivedSignal(configuration, simulationSteps + 1);
simulatedSignal = applyChannelIR(simulatedSignalRaw, channelTaps);
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
        z_hat_k = measurementFunction(x_k_k_1, configuration) / samplesTotal;
        
        if plotMeasures
            % ---- Plot routine -----
            plot(correlatorTaps, real(z_k));
            hold on;
            plot(correlatorTaps, real(z_hat_k));
            plot(correlatorTaps, imag(z_k));
            plot(correlatorTaps, imag(z_hat_k));
            hold off;
            ylabel('Real and Imag parts of z_k and z_k_hat');
            xlabel('Correlator tap');
            legend({'Real $z[k]$', 'Real $\hat{z}[k]$', 'Imag $z[k]$', 'Imag $\hat{z}[k]$'}, 'Interpreter','latex');
            pause(0.1);
        end
        
        % Compute Jacobian
        % Delay term now follows Φ_pp(ετ + (l - m)Ts) as in the analytical model.
        delayJacobian = delayJacobianFunctionSimplified( ...
            x_k_k_1(1), ...
            x_k_k_1(5:end), ...
            1 / configuration.samplingFrequency, ...
            q, ...
            1 / configuration.chippingFrequency, ...
            exp(1j * x_k_k_1(2))  ...
        );
        phaseJacobian = 1j * z_hat_k;
        dopplerJacobian = zeros(2*q + 1, 2);
        channelOrder = numel(x_k_k_1(5:end)) - 1;
        channelWeightsJacobian = exp(1j * x_k_k_1(2)) .* ...
            getShiftedCorrelations(x_k_k_1(1), q, configuration, channelOrder) / samplesTotal;
        jacobian = [delayJacobian ...
            phaseJacobian ...
            dopplerJacobian ...
            channelWeightsJacobian];
        
        % Compute Kalman Gain
        K_k = P_k_k_1 * jacobian'...
            *((jacobian * P_k_k_1 * jacobian' + R) \ eye(C));
        kalmanGainRecord(:,:, k) = K_k;

        % Obtain the innovation
        innovation = z_k - z_hat_k;
        innovationRecord(:, k) = innovation;

        % EKF's state update
        x_k_k = x_k_k_1 + K_k * innovation;
        x_k_k(WienerStatesSelection) = real(x_k_k(WienerStatesSelection));
        
        % EKF's covariance matrix update
        P_k_k = (eye(q + 1 + 4) - K_k*jacobian) * P_k_k_1;
        
        % LQG control vector computation
        u_LQG = L * x_k_k(WienerStatesSelection);
    else
        % Initialization procedure
        x_k_k = x_k_k_1;
        P_k_k = P_k_k_1;
    end

    % EKF's Projection Ahead Step
    x_k_k_1 = F * x_k_k;  
    x_k_k_1(WienerStatesSelection) = real(x_k_k_1(WienerStatesSelection));
    P_k_k_1 = F * P_k_k * F' + Q;

    errorStateRecord(:, k) = x_k_k_1(1:4);
    channelStateRecord(:, k) = x_k_k_1(5:end);
end
%% Plots

lineWidth = 2;
fontSize = 13;

epochVector = 1:simulationSteps;

% Observe the STD of the innovation sequence time series
innovationStdRecord = std(innovationRecord,1,1);
figure(Name="STD of the innovations", NumberTitle="off");
plot(epochVector, innovationStdRecord, 'LineWidth', lineWidth);
ylabel("Standard deviation of the innovations");
xlabel("Epochs (Simulation Steps)");

% Observe the innovation sequence time series
figure(Name="Middle tap of the innovation sequence", NumberTitle="off");
hold on;
plot(epochVector, real(innovationRecord(middleSample,:)), 'LineWidth', lineWidth);
plot(epochVector, imag(innovationRecord(middleSample,:)), 'LineWidth', lineWidth);
legend({"Real", "Imaginary"});
ylabel("Innovation sequence of the middle tap");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="Delay Estimation", NumberTitle="off");
hold on;
plot(epochVector, LQGStateRecord(1,:), 'LineWidth', lineWidth);
plot(epochVector, LOSDelay(epochVector*4096), 'LineWidth', lineWidth);
legend({"LQG's estimated delay", "True delay"});
ylabel("Delay estimate");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="Delay Error State", NumberTitle="off");
hold on;
plot(errorStateRecord(1, :), 'LineWidth', lineWidth);
plot(epochVector, zeros(1, simulationSteps), 'LineWidth', lineWidth);
legend({"EKF's estimated delay error", "Zero line"});
ylabel("Delay error estimate");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="Phase Estimation", NumberTitle="off");
hold on;
plot(epochVector, LQGStateRecord(2,:), 'LineWidth', lineWidth);
plot(epochVector, LOSPhase(epochVector*4096), 'LineWidth', lineWidth);
legend({"LQG's estimated phase", "True Phase"});
ylabel("Phase estimate");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="Phase Error State", NumberTitle="off");
hold on;
plot(errorStateRecord(2, :), 'LineWidth', lineWidth);
plot(epochVector, zeros(1, simulationSteps), 'LineWidth', lineWidth);
legend({"EKF's estimated phase error", "Zero line"});
ylabel("Phase error estimate");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="Doppler Estimation", NumberTitle="off");
hold on;
plot(epochVector, LQGStateRecord(3,:), 'LineWidth', lineWidth);
plot(epochVector, 2*pi*configuration.dopplerProfile(2) * ones(1,length(epochVector)), 'LineWidth', lineWidth);
legend({"LQG's estimated Doppler frequency", "True Doppler frequency"});
ylabel("Doppler estimate");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="Doppler Error State", NumberTitle="off");
hold on;
plot(errorStateRecord(3, :), 'LineWidth', lineWidth);
plot(epochVector, zeros(1, simulationSteps), 'LineWidth', lineWidth);
ylabel("Doppler error estimate");
xlabel("Epochs (Simulation Steps)");
legend({"EKF's estimated Doppler error", "Zero line"});
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

% Channel taps — boxplots (Real) with true (scatter)
figure(Name="Channel Taps (Real) — Boxplot vs True", NumberTitle="off"); hold on;
boxplot(real(channelStateRecord.'));
m = min(size(channelStateRecord,1), numel(channelTaps));
scatter(1:m, real(channelTaps(1:m)), 36, 'b', 'filled', 'DisplayName','True');
xticklabels(arrayfun(@num2str, 0:size(channelStateRecord,1)-1, 'UniformOutput', false));
xlabel("Tap index (samples)"); ylabel("Real part");
set(gca, "FontSize", fontSize); grid on; legend('Location','best'); hold off;

% Channel taps — boxplots (Imag) with true (scatter)
figure(Name="Channel Taps (Imag) — Boxplot vs True", NumberTitle="off"); hold on;
boxplot(imag(channelStateRecord.')); 
scatter(1:m, imag(channelTaps(1:m)), 36, 'b', 'filled', 'DisplayName','True');
xticklabels(arrayfun(@num2str, 0:size(channelStateRecord,1)-1, 'UniformOutput', false));
xlabel("Tap index (samples)"); ylabel("Imag part");
set(gca, "FontSize", fontSize); grid on; legend('Location','best'); hold off;
