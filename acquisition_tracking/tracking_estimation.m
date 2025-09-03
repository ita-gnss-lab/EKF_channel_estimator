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

varianceSquared = [1e-3 1e-3 1e-2 1e-2 1e-2];
stateTransitionCovariance = getCovarianceMatrix(...
    varianceSquared, ...
    epoch, ...
    configuration.carrierFrequency, ...
    numberOfTaps);

%% State History Vectors
delay = zeros(1, simulationSteps);
error = zeros(1, simulationSteps);

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
ECostMatrix = 0.1 * blkdiag(beta, 1, 1/epoch, 2/epoch^2);
UCostMatrix = 0.1 * blkdiag(beta, 1, 1/epoch, 2/epoch^2);

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
stateAPosteriori(1) = 0;
stateAPosteriori(5) = 1;

stateCovarianceMatrixAPosteriori = blkdiag(1e-6, (2*pi)^2/12, (50)^2/12, 0.2^2/12, 0.01 * eye(1 + numberOfTaps));

carrierState = [1.1e-4 configuration.dopplerProfile].';

% (Rodrigo: Compute the controlInput as L * stateAPosteriori)
controlInput = controlMatrix * stateAPosteriori(1:4);

%% Simulation
for k = 1 : simulationSteps
    %% Forward Step
    stateAPriori = stateTransitionMatrix * stateAPosteriori;  
    stateAPriori(carrierError) = real(stateAPriori(carrierError));
    stateCovarianceMatrixAPriori = stateTransitionMatrix * ...
        stateCovarianceMatrixAPosteriori * stateTransitionMatrix' ...
        + stateTransitionCovariance;

    error(k) = stateAPosteriori(1);

    %% Simulate Signal
    % (Rodrigo): Put this out of the loop
    [receivedSignal, time] = gnss_received_signal(configuration, epoch);
    samplesTotal = length(time);
    
    %% Carrier Removal
    
    % Update State 
    carrierState = ...
        carrierStateTransitionMatrix * carrierState + ...
        carrierCouplingMatrix * controlInput;

    delay(k) = carrierState(1);
    
    % Carrier Wipe-Off
    [totalPhaseAPriori, ~, ~] = get_LOS_dynamics(...
        time, ...
        carrierState(2:4).', ...
        configuration.carrierFrequency);
    carrierCorrection = exp(1j * totalPhaseAPriori);
    
    wipedSignal = receivedSignal .* conj(carrierCorrection);
    
    %% Multi-Correlator 
    delayAPriori = carrierState(1);
    delaysVector = delayAPriori + ...
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
plot(delay);
hold on;
plot(1e-4 * ones(1, simulationSteps));

figure(Name="Error Estimation", NumberTitle="off");
plot(error);
hold on;
plot(zeros(1, simulationSteps));


