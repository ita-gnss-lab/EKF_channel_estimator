clear;
load config_no_doppler.mat

%% Parameters
simulationSteps = 500;
numberOfTaps = 2;
totalChips = 1023;
epoch = totalChips / configuration.chippingFrequency;

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

%% Initial State
stateAPosteriori = zeros(numberOfTaps + 5, 1);
stateAPosteriori(1) = 2e-4;
stateAPosteriori(5) = 1;
stateCovarianceMatrixAPosteriori = blkdiag(1e-4, 100, 1e-5, 1e-5, eye(1 + numberOfTaps));

carrierState = [1e-3 configuration.dopplerProfile].';

controlInput = zeros(4, 1);

%% Simulation
for k = 1 : simulationSteps
    %% Forward Step
    stateAPriori = stateTransitionMatrix * stateAPosteriori;  
    stateAPriori(1) = real(stateAPriori(1));
    disp(stateAPriori(1));
    stateCovarianceMatrixAPriori = stateTransitionMatrix * ...
        stateCovarianceMatrixAPosteriori * stateTransitionMatrix.' ...
        + stateTransitionCovariance;

    %% Simulate Signal
    [receivedSignal, time] = gnss_received_signal(configuration, epoch);
    samplesTotal = length(time);
    
    %% Carrier Removal
    
    % Update State 
    carrierState = ...
        carrierStateTransitionMatrix * carrierState + ...
        carrierCouplingMatrix * controlInput;
    
    % Carrier Wipe-Off
    [totalPhaseAPriori, ~, ~] = get_LOS_dynamics(...
        time, ...
        stateAPriori(2:4).', ...
        configuration.carrierFrequency);
    carrierCorrection = exp(totalPhaseAPriori);
    
    wipedSignal = receivedSignal .* conj(carrierCorrection);
    
    %% Multi-Correlator 
    delayAPriori = stateAPriori(1);
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
    channelWeights = stateAPriori(5:end);
    channelWeightsJacobian = exp(stateAPriori(2)) * channelWeights.' ...
        .* getShiftedCorrelations(stateAPriori(1), numberOfTaps, configuration);

    jacobian = [delayJacobian ...
        phaseJacobian ...
        dopplerJacobian ...
        channelWeightsJacobian];
    
    %% Kalman Filter Estimation  
    
    kalmanGain = stateCovarianceMatrixAPriori * jacobian'...
        * inv(jacobian * stateCovarianceMatrixAPriori * jacobian' + noiseCovarianceMatrix);
    stateAPosteriori = stateAPriori + kalmanGain * (measurement - measurementEstimative);
    stateAPosteriori(1) = real(stateAPosteriori(1));
    stateCovarianceMatrixAPosteriori = (eye() - kalmanGain*jacobian) * ...
        stateCovarianceMatrixAPriori;
    
    [stablizedCostMatrix, controlMatrix, ~] = idare(carrierStateTransitionMatrix, ...
        carrierCouplingMatrix, ...
        ECostMatrix, ...
        UCostMatrix, ...
        [], []);
    
    errorStateAPosteriori = stateAPosteriori(1:4); 
    controlInput = controlMatrix * errorStateAPosteriori;

end