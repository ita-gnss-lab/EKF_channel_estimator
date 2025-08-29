load config_no_doppler.mat

%% Parameters
simulationSteps = 500;

%% Covariances 
% Convert CN0 from dB-Hz to linear scale
carrierToNoiseRatio = 10^(configuration.carrierToNoiseDensityRatio / 10);

% Compute the noise variance
termalNoiseVarianceSquared = configuration.samplingFrequency / carrierToNoiseRatio;

%% Simulation
for k = 1 : simulationSteps
    %% Forward Step
    stateAPriori = stateTransitionMatrix * stateAPosteriori;  
    stateCovarianceMatrixAPriori = stateTransitionMatrix * ...
        stateCovarianceMatrixAPosteriori * stateTransitionMatrix.' ...
        + Q;
    
    %% Simulate Signal
    
    [receivedSignal, time] = gnss_received_signal(configuration, epoch);
    
    %% Carrier Removal
    
    % Update State 
    carrierState = ...
        carrierStateTransitonMatrix * carrierState + ...
        carrierCouplingMatrix * controlInput;
    
    % Carrier Wipe-Off
    carrierCorrection = exp(totalPhaseAPriori);
    
    wipedSignal = receivedSignal * carrierCorrection';
    
    %% Multi-Correlator 
    
    delaysVector = delayAPriori + tapSpacing * ...
        (-numberOfTaps : 1 : numberOfTaps);
    correlatorBank = zeros(length(delaysVector), ...
        configuration.samplingFrequency * deltaTime);
    for i = 1:length(delaysVector)
        correlatorBank(i,:) = reference_signal(configuration.satellite, ...
                    delaysVector(i), ...
                    configuration.chippingFrequency, ...
                    configuration.sampplingFrequency, ...
                    1)';
    end
    
    measurement = correlatorBank * wipedSignal / samplesTotal;
    measurementEstimative = measurementFunction(stateAPriori);

    noiseCovarianceMatrix = ...
        (termalNoiseVarianceSquared / samplesTotal.^2) * ...
        (correlatorBank * correlatorBank.');
    
    %% Compute Jacobian
    delayJacobian = delayJacobianFunction(stateAPriori);
    phaseJacobian = 1j * measurementEstimative;
    dopplerJacobian = zeros(2 * numberOfTaps + 1, 2);
    channelWeightsJacobian = exp(totalPhaseError) * channelWeights ...
        .* getShiftedCorrelations(delayError, numberOfTaps, samplingPeriod);
    
    %% Kalman Filter Estimation  
    
    kalmanGain = stateCovarianceMatrixAPriori * jacobian' \ ...
        (jacobian * stateCovarianceMatrixAPriori * jacobian' + noiseCovarianceMatrix);
    stateAPosteriori = stateAPriori + kalmanGain * (measurement - measurementEstimative);
    stateCovarianceMatrixAPosteriori = (eye() - kalmanGain*jacobian) * ...
        stateCovarianceMatrixAPriori;
    
    [stablizedCostMatrix, controlMatrix, ~] = idare(transitionMatrix, ...
        couplingMatrix, ...
        ECostMatrix, ...
        UCostMatrix, ...
        [], []);
    
    errorStateAPosteriori = stateAPosteriori(selection); 
    controlInput = controlMatrix * errorStateAPosteriori;

end