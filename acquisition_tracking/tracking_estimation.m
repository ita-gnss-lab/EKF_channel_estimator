clearvars; clc; close all;

addpath(genpath(fullfile("..", "..","EKF_channel_estimator")));

load config_no_doppler.mat
rng(26437226);

%% Parameters
simulationSteps = 500;
q = 2;
C = 2*q + 1;
middleSample = q + 1;
epoch = configuration.totalChips / configuration.chippingFrequency;

% NOTE: These were not being used
WienerStatesSelection = 1:4;

%% Covariances 
% Convert CN0 from dB-Hz to linear scale
carrierToNoiseRatioLinear = 10^(configuration.carrierToNoiseDensityRatio / 10);
% Compute the noise variance
thermalNoiseVarianceSquared = configuration.samplingFrequency / carrierToNoiseRatioLinear;

sigma2Vec = [1e-4 1e-8 1e-6 1e-6 0];
Q = getStateCovarianceMatrix(...
    sigma2Vec, ...
    epoch, ...
    configuration.carrierFrequency, ...
    q);

%% State History Vectors
LQGStateRecord = zeros(4, simulationSteps);
errorStateRecord = zeros(4, simulationSteps);
channelStateRecord = zeros(3, simulationSteps);
innovationRecord = zeros(C, simulationSteps);

%% Transition Matrices
[F_W, F_H] = getModelTransitionMatrix( epoch, ...
    configuration.carrierFrequency, q);

F = blkdiag(F_W, F_H);

%% Cost Functions
beta = 1 / (2 * pi * configuration.carrierFrequency);
relation = 0.5;
T_e =  relation * blkdiag(beta, 1, 1/epoch, 2/epoch^2);
T_u =  blkdiag(beta, 1, 1/epoch, 2/epoch^2);

%% Coupling Matrix for Control Signal
B_LQG = eye(4);

%% IDARE Solution
[~, L, ~] = idare(F_W, B_LQG, T_e, T_u, [], []);

%% Initial State
x_hat_k_k = zeros(q + 5, 1);
x_hat_k_k(5) = 1;

% HACK: I zeroed this initial covariance matrix to my analysis about the
% phase estimation.
channelCovarianceMatrix = 0 * eye(1 + q); %0.000001 * eye(1 + q);
channelCovarianceMatrix(1,1) = 0; % 0.001;  

P_k_k = blkdiag(1e-9, (2*pi)^2/12, 0.0001*(50)^2/12, 0, channelCovarianceMatrix); 
% P_k_k = blkdiag(0, 0, 0, 0, zeros(1 + q));

phaseError = 0;
x_LQG_k = [1.01e-4, ...
    configuration.dopplerProfile(1) + phaseError, ...
    2*pi*configuration.dopplerProfile(2:end)].';

u_LQG = L * x_hat_k_k(WienerStatesSelection);

%% Simulate Signal
configuration.addNoise = false;
[simulatedSignal, ~, LOSPhase, LOSDelay] = gnssReceivedSignal(configuration, simulationSteps + 1);
samplesTotal = epoch*configuration.samplingFrequency;

%% Simulation
plotMeasures = false;
for k = 1 : simulationSteps
    %% Forward Step
    x_k_k_1 = F * x_hat_k_k;  
    x_k_k_1(WienerStatesSelection) = real(x_k_k_1(WienerStatesSelection));
    P_k_k_1 = F * ...
        P_k_k * F' ...
        + Q;

    errorStateRecord(:, k) = x_hat_k_k(1:4);

    %% Signal 
    receivedSignal = simulatedSignal(((k - 1) * samplesTotal + 1: k * samplesTotal));
    time = 1 / configuration.samplingFrequency : 1 / configuration.samplingFrequency : epoch;

    %% Carrier Removal
    
    % Update State 
    x_LQG_k = F_W * x_LQG_k - B_LQG * u_LQG;

    LQGStateRecord(:, k) = x_LQG_k(1:4);
    
    % Carrier Wipe-Off
    [totalPhaseAPriori, actualUsedDelay, ~] = get_LOS_dynamics(...
        time, ...
         x_LQG_k(2:4).', ...  % configuration.dopplerProfile
        configuration.carrierFrequency);
    carrierCorrection = exp(1j * totalPhaseAPriori);
    
    wipedSignal = receivedSignal .* conj(carrierCorrection);
    
    %% Multi-Correlator 

    % HACK(Rodrigo): I used the (known and fixed) value of 1e-4 to debug
    % buildCorrelatorBank.
    delayAPriori = x_LQG_k(1); %1e-4 

    % NOTE(Rodrigo): This buildCorrelatorBank is a new function that i 
    % created to build the correlator bank. In my understanding, we were
    % building the correlator bank in an incorrect manner previously, since
    % the samples delay were being computed individually for each code
    % replica. I think it is better to compute a single delay in samples
    % and then use it to build the correlator bank by shifting a code
    % replica using integer values of samples delay.
    correlatorBank = buildCorrelatorBank(configuration, delayAPriori, q);

    % delaysVector = delayAPriori + ...
    %     1 / (configuration.chippingFrequency * q) * ...
    %     (-q : 1 : q);
    % correlatorBank = zeros(length(delaysVector), ...
    %     samplesTotal);
    % for i = 1:length(delaysVector)
    %     correlatorBank(i,:) = reference_signal(configuration, ...
    %                 delaysVector(i))';
    % end
    
    % NOTE(Rodrigo): Put now a debug in measurementEstimate and plot
    % z_k. You can now see a perfect triangle, as we would expect.
    z_k = correlatorBank * wipedSignal / samplesTotal;

    % HACK(Rodrigo): I'm aritfically inputing the perfect version of, 
    % StateAPriori, so we can further modify the measurementFunction
    % function so that the shape of the correlation matches what we are
    % getting from the buildCorrelatorBank function.
    % z_hat_k = measurementFunction([zeros(4,1);1;zeros(q, 1)], ...
    %     configuration) / samplesTotal;
    z_hat_k = measurementFunction(x_k_k_1, ...
        configuration) / samplesTotal;
    
    if plotMeasures
        % ---- Plot routine -----
        plot(abs(z_k));
        hold on;
        plot(abs(z_hat_k));
        hold off;
        pause(0.1)
    end

    R = (thermalNoiseVarianceSquared / samplesTotal.^2) * ...
        (correlatorBank * correlatorBank.');
    
    %% Compute Jacobian
    % delayJacobian = delayJacobianFunction(x_k_k_1, configuration) / samplesTotal;
    % delayJacobian = [1/configuration.chippingFrequency ; 1/configuration.chippingFrequency ; 0 ; -1/configuration.chippingFrequency; -1/configuration.chippingFrequency];
    delayJacobian = delayJacobianFunctionSimplified( ...
        x_k_k_1(1), ...
        x_k_k_1(5:end), ...
        1 / configuration.samplingFrequency, ...
        q, ...
        1 / configuration.chippingFrequency, ...
        sqrt(1) * exp(1j * x_k_k_1(2))  ...
    );
    phaseJacobian = 1j * z_hat_k;
    dopplerJacobian = zeros(2*q + 1, 2);
    channelWeightsJacobian = exp(1j * x_k_k_1(2)) .* ...
        getShiftedCorrelations(x_k_k_1(1), q, configuration) / samplesTotal;

    jacobian = [delayJacobian ...
        phaseJacobian ...
        dopplerJacobian ...
        channelWeightsJacobian];
    
    %% Kalman Filter Estimation  
    
    % HACK: The inversion of the matrix is hard-coded to use an eye(5) (and
    % eye(7) somewhere else) matrix. We need to make it more general later.
    K_k = P_k_k_1 * jacobian'...
        *((jacobian * P_k_k_1 * jacobian' + R) \ eye(C));

    innovation = z_k - z_hat_k;
    innovationRecord(:, k) = innovation;
    x_hat_k_k = x_k_k_1 + K_k * innovation;
    x_hat_k_k(WienerStatesSelection) = real(x_hat_k_k(WienerStatesSelection));
    % x_hat_k_k(5:end) = [1 0 0].';
    P_k_k = (eye(q + 1 + 4) - K_k*jacobian) * ...
        P_k_k_1;
    
    %% Control Signal Computation
    
    u_LQG = L * x_hat_k_k(WienerStatesSelection);
end
%% Plots

epochVector = 1:simulationSteps;

% Observe the STD of the innovation sequence time series
innovationStdRecord = std(innovationRecord,1,1);
figure(Name="STD of the innovations", NumberTitle="off");
plot(epochVector, innovationStdRecord);
ylabel("Standard deviation of the innovations");
xlabel("Epochs (Simulation Steps)");

% Observe the innovation sequence time series
figure(Name="Middle tap of the innovation sequence", NumberTitle="off");
plot(epochVector, abs(innovationRecord(middleSample,:)));
ylabel("Innovation sequence of the middle tap");
xlabel("Epochs (Simulation Steps)");

figure(Name="Delay Estimation", NumberTitle="off");
plot(epochVector, LQGStateRecord(1,:));
hold on;
plot(epochVector, LOSDelay(epochVector*4000));
ylabel("Delay estimate");
xlabel("Epochs (Simulation Steps)");

figure(Name="Delay Error State", NumberTitle="off");
plot(errorStateRecord(1, :));
hold on;
plot(epochVector, zeros(1, simulationSteps));
ylabel("Delay error estimate");
xlabel("Epochs (Simulation Steps)");

figure(Name="Phase Estimation", NumberTitle="off");
plot(epochVector, LQGStateRecord(2,:));
hold on;
plot(epochVector, LOSPhase(400*epochVector));
ylabel("Doppler estimate");
xlabel("Epochs (Simulation Steps)");

figure(Name="Phase Error State", NumberTitle="off");
plot(errorStateRecord(2, :));
hold on;
plot(epochVector, zeros(1, simulationSteps));
ylabel("Doppler error estimate");
xlabel("Epochs (Simulation Steps)");

figure(Name="Doppler Estimation", NumberTitle="off");
plot(epochVector, LQGStateRecord(3,:));
hold on;
plot(epochVector, 2*pi*configuration.dopplerProfile(2) * ones(1,length(epochVector)));
ylabel("Doppler estimate");
xlabel("Epochs (Simulation Steps)");

figure(Name="Doppler Error State", NumberTitle="off");
plot(errorStateRecord(3, :));
hold on;
plot(epochVector, zeros(1, simulationSteps));
ylabel("Doppler error estimate");
xlabel("Epochs (Simulation Steps)");
