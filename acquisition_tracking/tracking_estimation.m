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

sigma2WVec = [1e-4 1e-6 1e-6 1e-6 1e-6];
Q = getStateCovarianceMatrix(...
    sigma2WVec, ...
    epoch, ...
    configuration.carrierFrequency, ...
    q);

%% State History Vectors
LQGStateRecord = zeros(4, simulationSteps);
errorStateRecord = zeros(4, simulationSteps);
channelStateRecord = zeros(3, simulationSteps);
innovationRecord = zeros(C, simulationSteps);

%% Transition Matrices
[F_W, F_H] = ...
    getModelTransitionMatrix(...
    epoch, ...
    configuration.carrierFrequency, ...
    q);

F = blkdiag(...
    F_W,...
    F_H);

%% Cost Functions
beta = 1 / (2 * pi * configuration.carrierFrequency);
relation = 0.3;
T_e =  relation * blkdiag(beta, 1, 1/epoch, 2/epoch^2);
T_u =  blkdiag(beta, 1, 1/epoch, 2/epoch^2);

%% Coupling Matrix for Control Signal
B_LQG = eye(4);

%% IDARE Solution
[~, L, ~] = idare(F_W, ...
    B_LQG, ...
    T_e, ...
    T_u, ...
    [], []);

%% Initial State
x_hat_k_k = zeros(q + 5, 1);
x_hat_k_k(5) = 1;

channelCovarianceMatrix = 0.000001 * eye(1 + q);
channelCovarianceMatrix(1,1) = 0.01;  

P_k_k = blkdiag(1e-8, (2*pi)^2/12, (50)^2/12, 0.2^2/12, channelCovarianceMatrix);
% stateCovarianceMatrixAPosteriori = blkdiag(0, 0, 0, 0, zeros(1 + numberOfTaps));

x_LQG_k = [0.995e-4 configuration.dopplerProfile].';

u_LQG = L * x_hat_k_k(WienerStatesSelection);

%% Simulate Signal
configuration.addNoise = true;
[simulatedSignal, ~] = gnssReceivedSignal(configuration, simulationSteps + 1);
samplesTotal = epoch*configuration.samplingFrequency;

%% Simulation
plotMeasures = false;
for k = 2 : simulationSteps
    %% Forward Step
    stateAPriori = F * x_hat_k_k;  
    stateAPriori(WienerStatesSelection) = real(stateAPriori(WienerStatesSelection));
    stateCovarianceMatrixAPriori = F * ...
        P_k_k * F' ...
        + Q;

    errorStateRecord(:, k) = x_hat_k_k(1:4);

    %% Signal 
    receivedSignal = simulatedSignal(((k - 1) * samplesTotal + 1: k * samplesTotal));
    time = 1 / configuration.samplingFrequency : 1 / configuration.samplingFrequency : epoch;

    %% Carrier Removal
    
    % Update State 
    x_LQG_k = ...
        F_W * x_LQG_k - ...
        B_LQG * u_LQG;

    LQGStateRecord(:, k) = x_LQG_k(1:4);
    
    % Carrier Wipe-Off
    [totalPhaseAPriori, actualUsedDelay, ~] = get_LOS_dynamics(...
        time, ...
        configuration.dopplerProfile, ...  %carrierState(2:4).'
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
    % measurement. You can now see a perfect triangle, as we would expect.
    measurement = correlatorBank * wipedSignal / samplesTotal;

    % HACK(Rodrigo): I'm aritfically inputing the perfect version of, 
    % StateAPriori, so we can further modify the measurementFunction
    % function so that the shape of the correlation matches what we are
    % getting from the buildCorrelatorBank function.
    % measurementEstimative = measurementFunction([zeros(4,1);1;zeros(q, 1)], ...
    %     configuration) / samplesTotal;
    measurementEstimative = measurementFunction(stateAPriori, ...
        configuration) / samplesTotal;
    
    if plotMeasures
        % ---- Plot routine -----
        plot(abs(measurement));
        hold on;
        plot(abs(measurementEstimative));
        hold off;
        pause(0.1)
    end

    noiseCovarianceMatrix = ...
        (thermalNoiseVarianceSquared / samplesTotal.^2) * ...
        (correlatorBank * correlatorBank.');
    
    %% Compute Jacobian
    % delayJacobian = delayJacobianFunction(stateAPriori, configuration) / samplesTotal;
    % delayJacobian = [1/configuration.chippingFrequency ; 1/configuration.chippingFrequency ; 0 ; -1/configuration.chippingFrequency; -1/configuration.chippingFrequency];
    delayJacobian = delayJacobianFunctionSimplified( ...
        stateAPriori(1), ...
        stateAPriori(5:end), ...
        1 / configuration.samplingFrequency, ...
        q, ...
        1 / configuration.chippingFrequency, ...
        sqrt(1) * exp(1j * stateAPriori(2))  ...
    );
    phaseJacobian = 1j * measurementEstimative;
    dopplerJacobian = zeros(2*q + 1, 2);
    channelWeightsJacobian = exp(1j * stateAPriori(2)) .* ...
        getShiftedCorrelations(stateAPriori(1), q, configuration) / samplesTotal;

    jacobian = [delayJacobian ...
        phaseJacobian ...
        dopplerJacobian ...
        channelWeightsJacobian];
    
    %% Kalman Filter Estimation  
    
    % HACK: The inversion of the matrix is hard-coded to use an eye(5) (and
    % eye(7) somewhere else) matrix. We need to make it more general later.
    kalmanGain = stateCovarianceMatrixAPriori * jacobian'...
        *((jacobian * stateCovarianceMatrixAPriori * jacobian' + noiseCovarianceMatrix) \ eye(C));

    innovation = measurement - measurementEstimative;
    innovationRecord(:, k) = innovation;
    x_hat_k_k = stateAPriori + kalmanGain * innovation;
    x_hat_k_k(WienerStatesSelection) = real(x_hat_k_k(WienerStatesSelection));
    % x_hat_k_k(5:end) = [1 0 0].';
    P_k_k = (eye(q + 1 + 4) - kalmanGain*jacobian) * ...
        stateCovarianceMatrixAPriori;
    
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
plot(1e-4 * ones(1, simulationSteps));
ylabel("Delay estimate");
xlabel("Epochs (Simulation Steps)");

figure(Name="Delay Error State", NumberTitle="off");
plot(errorStateRecord(1, :));
hold on;
plot(epochVector, zeros(1, simulationSteps));
ylabel("Delay error estimate");
xlabel("Epochs (Simulation Steps)");
