clear;
load config_no_doppler.mat
rng(26437226);

%% Parameters
simulationSteps = 500;
numberOfTaps = 2;
totalChips = 1023;
epoch = totalChips / configuration.chippingFrequency;
carrierError = 1:4;
channelWeights = 5 + 0:numberOfTaps;

%% Covariances 
% Convert CN0 from dB-Hz to linear scale
carrierToNoiseRatio = 10^(configuration.carrierToNoiseDensityRatio / 10);
% Compute the noise variance
termalNoiseVarianceSquared = configuration.samplingFrequency / carrierToNoiseRatio;

varianceSquared = [1e-2 1e-3 1e-2 1e-2 1e-2];
stateTransitionCovariance = getCovarianceMatrix(...
    varianceSquared, ...
    epoch, ...
    configuration.carrierFrequency, ...
    numberOfTaps);

%% State History Vectors
carrierStateRecord = zeros(4, simulationSteps);
errorStateRecord = zeros(4, simulationSteps);

%% Transition Matrices
[carrierStateTransitionMatrix, channelStateTransitionMatrix] = ...
    getModelTransitionMatrix(...
    epoch, ...
    configuration.carrierFrequency, ...
    numberOfTaps);

stateTransitionMatrix = blkdiag(...
    carrierStateTransitionMatrix,...
    channelStateTransitionMatrix);

%% Cost Functions
beta = 1 / (2 * pi * configuration.carrierFrequency);
relation = 1.1;
ECostMatrix =  relation * blkdiag(beta, 1, 1/epoch, 2/epoch^2);
UCostMatrix =  blkdiag(beta, 1, 1/epoch, 2/epoch^2);

%% Coupling Matrix for Control Signal
carrierCouplingMatrix = eye(4);

%% IDARE Solution
[stablizedCostMatrix, controlMatrix, ~] = idare(carrierStateTransitionMatrix, ...
    carrierCouplingMatrix, ...
    ECostMatrix, ...
    UCostMatrix, ...
    [], []);

%% Initial State
stateAPosteriori = zeros(numberOfTaps + 5, 1);
stateAPosteriori(1) = -1e-6;
stateAPosteriori(5) = 1;

stateCovarianceMatrixAPosteriori = blkdiag(1e-6, (2*pi)^2/12, (50)^2/12, 0.2^2/12, 0.01 * eye(1 + numberOfTaps));

carrierState = [1.1e-4 configuration.dopplerProfile].';

controlInput = controlMatrix * stateAPosteriori(carrierError);

%% Simulate Signal
% (Rodrigo): Put this out of the loop
[simulatedSignal, totalTime] = gnss_received_signal(configuration, epoch*(simulationSteps + 1));
samplesTotal = epoch*configuration.samplingFrequency + 1;

%% Simulation
for k = 2 : simulationSteps
    %% Forward Step
    stateAPriori = stateTransitionMatrix * stateAPosteriori;  
    stateAPriori(carrierError) = real(stateAPriori(carrierError));
    stateCovarianceMatrixAPriori = stateTransitionMatrix * ...
        stateCovarianceMatrixAPosteriori * stateTransitionMatrix' ...
        + stateTransitionCovariance;

    errorStateRecord(:, k) = stateAPosteriori(1:4);

    %% Signal 
    receivedSignal = simulatedSignal(((k - 1) * samplesTotal + 1: k * samplesTotal));
    time = [0 : 1 / configuration.samplingFrequency : epoch];

    %% Carrier Removal
    
    % Update State 
    carrierState = ...
        carrierStateTransitionMatrix * carrierState + ...
        carrierCouplingMatrix * controlInput;
    disp('value');
    disp(carrierState(2));
    disp('error');
    disp(stateAPosteriori(2));

    carrierStateRecord(:, k) = carrierState(1:4);
    
    % Carrier Wipe-Off
    [totalPhaseAPriori, actualUsedDelay, ~] = get_LOS_dynamics(...
        time, ...
        carrierState(2:4).', ...
        configuration.carrierFrequency);
    carrierCorrection = exp(1j * totalPhaseAPriori);
    
    wipedSignal = receivedSignal .* conj(carrierCorrection);
    
    %% Multi-Correlator 
    delayAPriori = carrierState(1);
    delaysVector = actualUsedDelay(1) + ...
        1 / (configuration.chippingFrequency * numberOfTaps) * ...
        (-numberOfTaps : 1 : numberOfTaps);
    correlatorBank = zeros(length(delaysVector), ...
        samplesTotal);
    for i = 1:length(delaysVector)
        correlatorBank(i,:) = reference_signal(configuration.satellite, ...
                    delaysVector(i), ...
                    configuration.chippingFrequency, ...
                    configuration.samplingFrequency, ...
                    1)';
    end
    
    measurement = correlatorBank * wipedSignal / samplesTotal;
    measurementEstimative = measurementFunction(stateAPriori, ...
        configuration) / samplesTotal;
    plot(abs(measurement));
    %hold on;
    %plot(abs(measurementEstimative));
    %hold off;
    pause(0.1)

    noiseCovarianceMatrix = ...
        (termalNoiseVarianceSquared / samplesTotal.^2) * ...
        (correlatorBank * correlatorBank.');
    
    %% Compute Jacobian
    delayJacobian = delayJacobianFunction(stateAPriori, configuration);
    phaseJacobian = 1j * measurementEstimative;
    dopplerJacobian = zeros(2*numberOfTaps + 1, 2);
    channelWeightsJacobian = exp(1j * stateAPriori(2)) .* ...
        getShiftedCorrelations(stateAPriori(1), numberOfTaps, configuration);

    jacobian = [delayJacobian ...
        phaseJacobian ...
        dopplerJacobian ...
        channelWeightsJacobian];
    
    %% Kalman Filter Estimation  
    
    % HACK: The inversion of the matrix is hard-coded to use an eye(5) (and
    % eye(7) somewhere else) matrix. We need to make it more general later.
    kalmanGain = stateCovarianceMatrixAPriori * jacobian'...
        *((jacobian * stateCovarianceMatrixAPriori * jacobian' + noiseCovarianceMatrix) \ eye(5));
    stateAPosteriori = stateAPriori + kalmanGain * (measurement - measurementEstimative);
    stateAPosteriori(carrierError) = real(stateAPosteriori(carrierError));
    stateCovarianceMatrixAPosteriori = (eye(7) - kalmanGain*jacobian) * ...
        stateCovarianceMatrixAPriori;

    
    %% Control Signal Computation
    
    controlInput = controlMatrix * stateAPosteriori(carrierError);

end

figure(Name="Delay Estimation", NumberTitle="off");
plot(carrierStateRecord(1,:));
hold on;
plot(1e-4 * ones(1, simulationSteps));

figure(Name="Error Estimation", NumberTitle="off");
plot(errorStateRecord(1, :));
hold on;
plot(zeros(1, simulationSteps));


