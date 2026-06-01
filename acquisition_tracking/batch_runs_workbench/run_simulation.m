function run_simulation(experimentName, simulation_config, filter_config, seed)
    % RUN_EKF_BATCH Runs the EKF channel estimator simulation for batch processing.
    % Inputs:
    %   experimentName - String, prefix for the saved .mat files
    %   configuration  - Struct containing simulation parameters
    
    % Ensure path is set (you may want to handle this in your main batch script instead)
    addpath(genpath(fullfile("..", "..","EKF_channel_estimator")));
    
    rng(seed);
    samplesPerChip = simulation_config.samplesPerChip;
    simulation_config.samplingFrequency = samplesPerChip * simulation_config.chippingFrequency;
    
    %% Stress-test truth
    numberOfCausalTruthTaps = filter_config.numberOfCausalTruthTaps;
    truthTapOrder = 0:(numberOfCausalTruthTaps - 1);
    
    % If your batch caller doesn't define tdl_channel, define it here:
    if ~isfield(simulation_config, 'tdl_channel')
        simulation_config.tdl_channel = 0.8 * exp((-0.85 + 1j*pi/5) * truthTapOrder);
    end
    
    %% Parameters
    simulationSteps = 5000;
    useTapEnergyConstraint = false;
    useAcausalTaps = true;
    usePerfectFrozenTruthState = false;
    constraint_noise = 10^(-2.63);
    q = numberOfCausalTruthTaps;
    C = 2*q + 1;
    middleSample = q + 1;
    acausalTapStateIndices = 5:(4+q);
    acausalTapChannelIndices = 1:q;
    epoch = simulation_config.totalChips / simulation_config.chippingFrequency;
    simulation_config.correlatorHalfSpan = q;
    samplesTotal = epoch*simulation_config.samplingFrequency;
    timeSupport = (0:(samplesTotal - 1)).' * (1/simulation_config.samplingFrequency);
    beta = -1 / (2 * pi * simulation_config.carrierFrequency);
    WienerStatesSelection = 1:4;
    trueDelay = -simulation_config.dopplerProfile(1) / (2*pi*simulation_config.carrierFrequency);
    initialDelayErrorSamples = filter_config.initialDelayErrorSamples;
    initialDelayEstimate = trueDelay + initialDelayErrorSamples / simulation_config.samplingFrequency;
    initialDelayError = trueDelay - initialDelayEstimate;
    delayErrorStdSamples = max(2, abs(initialDelayErrorSamples));
    
    truthChannelAcausal = zeros(2*q + 1, 1);
    numberOfTruthTaps = min(numel(simulation_config.tdl_channel), q + 1);
    truthChannelAcausal(q + 1:q + numberOfTruthTaps) = simulation_config.tdl_channel(1:numberOfTruthTaps).';
    truthChannelState = truthChannelAcausal;
    
    %% Covariances 
    carrierToNoiseRatioLinear = 10^(simulation_config.carrierToNoiseDensityRatio / 10);
    thermalNoiseVarianceSquared = simulation_config.samplingFrequency / carrierToNoiseRatioLinear;
    sigma2Vec = [1e5 1e4 0 0 1e-2];
    Q = getStateCovarianceMatrix_acausal(sigma2Vec, epoch, beta, q);
    
    if ~useAcausalTaps
        Q(acausalTapStateIndices, :) = 0;
        Q(:, acausalTapStateIndices) = 0;
    end
    if usePerfectFrozenTruthState
        Q(:) = 0;
    end
    
    correlatorBank = buildCorrelatorBank(simulation_config, 0, q);
    R  = (thermalNoiseVarianceSquared / samplesTotal) * (correlatorBank * correlatorBank.');
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
    
    delayErrorStd = delayErrorStdSamples / simulation_config.samplingFrequency;
    phaseErrorStd = filter_config.phaseErrorStd;
    channelCovarianceMatrix = 1e-1 * eye(2*q + 1);
    channelCovarianceMatrix(q + 1, q + 1) = 1e-0;
    
    if usePerfectFrozenTruthState
        channelCovarianceMatrix(:) = 0;
    end
    if ~useAcausalTaps
        channelCovarianceMatrix(acausalTapChannelIndices, :) = 0;
        channelCovarianceMatrix(:, acausalTapChannelIndices) = 0;
    end
    
    P_k_k_1 = blkdiag(delayErrorStd^2, 10*phaseErrorStd^2, 0, 0, channelCovarianceMatrix);
    initialPhaseError = phaseErrorStd;
    x_LQG_k = [initialDelayEstimate, ...
        simulation_config.dopplerProfile(1) + initialPhaseError, ...
        2*pi*simulation_config.dopplerProfile(2), ...
        2*pi*simulation_config.dopplerProfile(3)].';
        
    u_LQG = L * x_k_k_1(WienerStatesSelection);
    u_LQG(3:4) = 0;
    
    trueTaps = simulation_config.tdl_channel;
    %% Simulate Signal
    [simulatedSignal, ~, LOSPhase, LOSDelay] = gnssReceivedSignal(simulation_config, simulationSteps + 1);
    increment = 0.1;
    %% Simulation
    for k = 1 : simulationSteps
        %% Signal 
        % Check if it is the first step
        if k == 1
            fprintf('Simulation started (Step 1/%d)\n', simulationSteps);
        % Check if k is a multiple of 10% of total steps
        elseif mod(k, floor(increment * simulationSteps)) == 0
            percentage = (k / simulationSteps) * 100;
            fprintf('Progress: %.0f%% complete (Step %d/%d)\n', percentage, k, simulationSteps);
        end
        receivedSignal = simulatedSignal(((k - 1) * samplesTotal + 1: k * samplesTotal));
        
        %% LQG Controller
        x_LQG_k = F_W * x_LQG_k + B_LQG * u_LQG;
        LQGStateRecord(:, k) = x_LQG_k(1:4);
        
        % Carrier Wipe-Off
        phi_T = x_LQG_k(2) + x_LQG_k(3) * timeSupport + 0.5 * x_LQG_k(4) * timeSupport.^2;
        d_k = exp(1j * phi_T);
        wipedSignal = receivedSignal .* conj(d_k);
        
        %% Kalman filter
        if k > 1
            % EKF's Update Step
            correlatorBank = buildCorrelatorBank(simulation_config, x_LQG_k(1), q);
            z_k = correlatorBank * wipedSignal / samplesTotal;
            constraint_value = sum(abs(x_k_k_1(other_taps)).^2)/((q-1)*abs(x_k_k_1(main_tap))^2);
            z_hat_k_aux = measurementFunction_acausal(x_k_k_1, simulation_config, q, x_LQG_k(1)) / samplesTotal;
            z_hat_k = z_hat_k_aux;
            
            if useTapEnergyConstraint
                z_k = [z_k; 0];
                z_hat_k = [z_hat_k; constraint_value];
            end
            constraintRecord(:, k) = constraint_value;
            
            % Compute Jacobian
            delayJacobian = delayJacobianFunctionSimplified_acausal( ...
                x_k_k_1(1), x_k_k_1(5:end), 1 / simulation_config.samplingFrequency, ...
                q, 1 / simulation_config.chippingFrequency, exp(1j * x_k_k_1(2)));
            phaseJacobian = 1j * z_hat_k_aux;
            dopplerJacobian = zeros(2*q + 1, 2);
            channelOrder = (numel(x_k_k_1(5:end)) - 1)/2;
            channelWeightsJacobian = exp(1j * x_k_k_1(2)) .* ...
                getShiftedCorrelations_acausal(x_k_k_1(1), q, simulation_config, channelOrder, x_LQG_k(1)) / samplesTotal;
            
            LOSParcel = -x_k_k_1(main_tap)*sum(abs(x_k_k_1(other_taps)).^2)/(abs(x_k_k_1(main_tap))^2);
            tapsParcel = (2/((q-1)*abs(x_k_k_1(main_tap))^2))*[x_k_k_1(4+1:4+q);LOSParcel;x_k_k_1(4+q+2:end)];
            constraintLine = [0 0 0 0 tapsParcel'];
            
            jacobian = [delayJacobian phaseJacobian dopplerJacobian channelWeightsJacobian];
            if useTapEnergyConstraint
                jacobian = [jacobian; constraintLine];
            end
            
            % Compute Kalman Gain
            K_k = P_k_k_1 * jacobian'*((jacobian * P_k_k_1 * jacobian' + R) \ eye(measurementDimension));
            kalmanGainRecord(:,:, k) = K_k;
            
            % Obtain the innovation
            innovation = z_k - z_hat_k;
            innovationRecord(:, k) = innovation;
            
            % EKF's state update
            x_k_k = x_k_k_1 + K_k * innovation;
            x_k_k(WienerStatesSelection) = real(x_k_k(WienerStatesSelection));
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
            u_LQG(3:4) = 0;
        else
            % Initialization procedure
            x_k_k = x_k_k_1;
            P_k_k = P_k_k_1;
        end
        
        % EKF's Projection Ahead Step
        x_k_k_1 = F * x_k_k;  
        x_k_k_1(WienerStatesSelection) = real(x_k_k_1(WienerStatesSelection));
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
    
    %% Compute Final Ground Truth Records (Math only, no plots)
    epochVector = 1:simulationSteps;
    truthSampleIndex = round(epochVector * samplesTotal);
    
    trueDelayRecord = LOSDelay(truthSampleIndex).';
    truePhaseRecord = LOSPhase(truthSampleIndex).';
    trueDopplerRecord = 2*pi*simulation_config.dopplerProfile(2) * ones(1, simulationSteps);
    
    trueDelayErrorRecord = trueDelayRecord - LQGStateRecord(1, :);
    truePhaseErrorRecord = truePhaseRecord - LQGStateRecord(2, :);
    trueDopplerErrorRecord = trueDopplerRecord - LQGStateRecord(3, :);
    
    innovationStdRecord = std(innovationRecord(1:C,:),1,1);
    
    VariationRecord = zeros(size(kalmanGainRecord, 1), simulationSteps);
    for i = 1:simulationSteps
        VariationRecord(:,i) = kalmanGainRecord(:,:,i)*innovationRecord(:,i);
    end

%% Save Record Variables
    % Define the folder path
    folderName = experimentName;
    
    % Check if the folder exists, if not, create it
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
    
    % Compile a list of variables you wish to extract
    recordVars = {
        'LQGStateRecord', 'errorStateRecord', 'channelStateRecord', ...
        'innovationRecord', 'kalmanGainRecord', 'constraintRecord', ...
        'VariationRecord', 'trueDelayRecord', 'truePhaseRecord', ...
        'trueDopplerRecord', 'trueDelayErrorRecord', ...
        'truePhaseErrorRecord', 'trueDopplerErrorRecord', 'innovationStdRecord', ...
        'trueTaps', 'simulation_config', 'filter_config', 'truthChannelState'
    };

    % Save each to an individual .mat file within the directory
    for i = 1:length(recordVars)
        varName = recordVars{i};
        % Construct the file path using the folderName
        fileName = fullfile(folderName, sprintf('%s_%s.mat', experimentName, varName));
        
        % Save the variable using the '-v7.3' flag
        save(fileName, varName, '-v7.3'); 
    end
end
