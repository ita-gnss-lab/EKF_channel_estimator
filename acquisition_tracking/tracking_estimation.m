load config_no_doppler.mat

%% Parameters
simulationSteps = 500;

%% Covariances 


%% Simulation
for k = 1 : simulationSteps
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
    stateCovarianceMatrixAPriori = (eye() - kalmanGain*jacobian) * ...
        stateCovarianceAPriori;
    
    [stablizedCostMatrix, controlMatrix, ~] = idare(transitionMatrix, ...
        couplingMatrix, ...
        ECostMatrix, ...
        UCostMatrix, ...
        [], []);
    
    errorStateAPosteriori = stateAPosteriori(selection); 
    controlInput = controlMatrix * errorStateAPosteriori;

end