function [delayRecord, delayTrue] = lqg_function(configuration, initialEstimate)
%LQG_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
    %% Parameters
    simulationSteps = configuration.steps;
    q = configuration.Ntaps;
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
    
    sigma2Vec = [1e-1 1e-1 1e-2 1e-3 0];
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
    channelStateRecord = zeros(3, simulationSteps);
    innovationRecord = zeros(C, simulationSteps);
    kalmanGainRecord = zeros(4 + q + 1,C, simulationSteps);
    
    %% Cost Function
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
    
    
    %% Initialization
    % NOTE: I changed from x_k_k to x_k_k_1, because, in fact the
    % initialization uses x[1|0].
    x_k_k_1 = zeros(q + 5, 1);
    x_k_k_1(5) = 1;
    
    % HACK: I zeroed this initial covariance matrix to my analysis about the
    % phase estimation.
    channelCovarianceMatrix = 0 * eye(1 + q); %0.000001 * eye(1 + q);
    channelCovarianceMatrix(1,1) = 0; % 0.001;  
    
    % NOTE: I changed from x_k_k to P_k_k_1, because, in fact the
    % initialization uses P[1|0].
    P_k_k_1 = blkdiag(1e-1, 1e-1, (50)^2/12, (0.1)^2/12, channelCovarianceMatrix); 
    % P_k_k_1 = blkdiag(1e-1, 0, 0, 0, zeros(1 + q));
    
    
    x_LQG_k = initialEstimate;
    % x_LQG_k = [1.005e-4, ...
    %     configuration.dopplerProfile(1) + phaseError, ...
    %     2*pi*(configuration.dopplerProfile(2) + DopplerError), ...
    %     2*pi*configuration.dopplerProfile(3)].';
    
    u_LQG = L * x_k_k_1(WienerStatesSelection);
    
    %% Simulate Signal
    configuration.addNoise = false;
    [simulatedSignal, ~, LOSPhase, LOSDelay] = gnssReceivedSignal(configuration, simulationSteps + 1);
    
    %% Simulation
    % NOTE: (Rodrigo): Changed the main loop to match algorithm 1 of my report.
    for k = 1 : simulationSteps
        disp(k*100/simulationSteps);
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
        errorStateRecord(:, k) = x_k_k(1:4);
    
    end

epochVector = 1:simulationSteps;

delayRecord = LQGStateRecord(1,:);
delayTrue = LOSDelay(epochVector*4096)';

end

