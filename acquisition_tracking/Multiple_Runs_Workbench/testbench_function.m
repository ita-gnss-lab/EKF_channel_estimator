function [innovationSum] = testbench_function(configuration, currentRun, terminationCondition, parentFolderName, initialCondition, covar_test)
    clc; close all;
    
    %addpath(genpath(fullfile("..", "..","EKF_channel_estimator")));
    if nargin < 3 || isempty(terminationCondition)
        terminationCondition = @(state) false; 
    end
    rng(26437226); 
    
    %% Parameters
    simulationSteps = 5000;
    constraint_noise = 10^(-4.44);%10.^[-2.63 -3.44 -4 -4.28 -4.49 -4.63 -4.85 -5.2 -5.37 -5.88];
    q = 5;
    C = 2*q + 1;
    middleSample = q + 1;
    epoch = configuration.totalChips / configuration.chippingFrequency;
    configuration.correlatorHalfSpan = q;
    samplesTotal = epoch*configuration.samplingFrequency;
    timeSupport = (0:(samplesTotal - 1)).' * (1/configuration.samplingFrequency);
    % NOTE: Verify if the signal should be + or - 
    beta = -1 / (2 * pi * configuration.carrierFrequency);
    WienerStatesSelection = 1:4;
    %simulationTaps = [1 0.1 0 0 0 0];
    
    %% Covariances 
    % Convert CN0 from dB-Hz to linear scale
    carrierToNoiseRatioLinear = 10^(configuration.carrierToNoiseDensityRatio / 10);
    % Compute the noise variancesamplesTotal.^2
    thermalNoiseVarianceSquared = configuration.samplingFrequency / carrierToNoiseRatioLinear;
    % [1e-1 1e-1 1e-2 1e-3 1e-4]
    sigma2Vec = [0 9.295e-5 7.877e-6 0 covar_test];
    Q = getStateCovarianceMatrix_acausal(...
        sigma2Vec, ...
        epoch, ...
        beta, ...
        q);
    correlatorBank = buildCorrelatorBank(configuration, 0, q);
    
    R  = (thermalNoiseVarianceSquared / samplesTotal) * ...
            (correlatorBank * correlatorBank.');
    R = [R zeros(2*q+1, 1); zeros(1, 2*q+1) constraint_noise]; 
    
    %% State History Vectors
    LQGStateRecord = zeros(4, simulationSteps);
    errorStateRecord = zeros(4, simulationSteps);
    channelStateRecord = zeros(2*q + 1, simulationSteps);
    innovationRecord = zeros(C+1, simulationSteps);
    kalmanGainRecord = zeros(4 + 2*q + 1,C+1, simulationSteps);
    constraintRecord = zeros(1, simulationSteps);
    PStateRecord = zeros(15, 15, simulationSteps);
    
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
    [F_W, F_H] = getModelTransitionMatrix_acausal(epoch, q, beta);
    
    %% Coupling Matrix for Control Signal
    B_LQG = eye(4);
    
    %% IDARE Solution
    [~, L, ~] = idare(F_W, B_LQG, T_e, T_u, [], []);
    
    %% Full transition matrix
    F = blkdiag(F_W, F_H);
    
    
    %% Initialization
    % NOTE: I changed from x_k_k to x_k_k_1, because, in fact the
    % initialization uses x[1|0].
    x_k_k_1 = zeros(4 + (2*q+1), 1);
    main_tap = false(size(x_k_k_1));
    main_tap(4 + q + 1) = true;
    x_k_k_1(main_tap) = initialCondition.mainTap;
    other_taps = true(size(x_k_k_1));
    other_taps(1:4) = false;
    other_taps(main_tap) = false;
    x_k_k_1(other_taps) = 0;
    
    % HACK: I zeroed this initial covariance matrix to my analysis about the
    % phase estimation.
    channelCovarianceMatrix = 1e-5 * eye(2*q + 1); %0.000001 * eye(1 + q);
    channelCovarianceMatrix(q + 1, q + 1) = 10*1e-5; % 0.001;  
  
    
    % NOTE: I changed from x_k_k to P_k_k_1, because, in fact the
    % initialization uses P[1|0].  1e-1, 0, (50)^2/12, (0.1)^2/12,
    P_k_k_1 = blkdiag((1/configuration.chippingFrequency)^2 / 3, pi^2/3, (100)^2 * 4 * pi^2 / 3, 0.2, channelCovarianceMatrix); 
    % P_k_k_1 = blkdiag(1e-1, 0, 0, 0, zeros(1 + q));
    
    
    phaseError = initialCondition.phaseError;
    DopplerError = initialCondition.DopplerError;
    delay = initialCondition.delay;
    x_LQG_k = [delay, ...
        configuration.dopplerProfile(1) + phaseError, ...
        2*pi*(configuration.dopplerProfile(2) + DopplerError), ...
        2*pi*configuration.dopplerProfile(3)].';
    
    u_LQG = L * x_k_k_1(WienerStatesSelection);
    
    %% Simulate Signal
    [simulatedSignal, ~, LOSPhase, LOSDelay] = gnssReceivedSignal(configuration, simulationSteps + 1);
    earlyBreak = 0;
    %simulatedSignal = applyChannelIR(simulatedSignal, simulationTaps);
    %% Simulation
    % NOTE: (Rodrigo): Changed the main loop to match algorithm 1 of my report.
    for k = 1 : simulationSteps
        %% Signal 
        receivedSignal = simulatedSignal(((k - 1) * samplesTotal + 1: k * samplesTotal));
        % NOTE (Thiago): corrects the samples that are not +-1
        receivedSignal = round(receivedSignal);
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
            % NOTE (Thiago): corrects the samples that are not +-1
            correlatorBank = round(correlatorBank);
            z_k = [correlatorBank * wipedSignal / samplesTotal; 0];
            constraint_value = sum(abs(x_k_k_1(other_taps)).^2)/((q-1)*abs(x_k_k_1(main_tap))^2);
            z_hat_k_aux = measurementFunction_acausal(x_k_k_1, configuration, q) / samplesTotal;
            z_hat_k = [z_hat_k_aux; constraint_value];
            constraintRecord(:, k) = constraint_value;
            
            
            % Compute Jacobian
            % Delay term now follows Φ_pp(ετ + (l - m)Ts) as in the analytical model.
            delayJacobian = delayJacobianFunctionSimplified_acausal( ...
                x_k_k_1(1), ...
                x_k_k_1(5:end), ...
                1 / configuration.samplingFrequency, ...
                q, ...
                1 / configuration.chippingFrequency, ...
                exp(1j * x_k_k_1(2)),  ...
                configuration.carrierFrequency...
            );
            phaseJacobian = 1j * z_hat_k_aux;
            dopplerJacobian = zeros(2*q + 1, 2);
            channelOrder = (numel(x_k_k_1(5:end)) - 1)/2;
            channelWeightsJacobian = exp(1j * x_k_k_1(2)) .* ...
                getShiftedCorrelations_acausal(x_k_k_1(1), q, configuration, channelOrder) / samplesTotal;
            % todo - add the jacobian of the constraint
            LOSParcel = -x_k_k_1(main_tap)*sum(abs(x_k_k_1(other_taps)).^2)/(abs(x_k_k_1(main_tap))^2);
            tapsParcel = (2/((q-1)*abs(x_k_k_1(main_tap))^2))*[x_k_k_1(4+1:4+q);LOSParcel;x_k_k_1(4+q+2:end)];
            constraintLine = [0 0 0 0 tapsParcel'];
            jacobian = [delayJacobian ...
                phaseJacobian ...
                dopplerJacobian ...
                channelWeightsJacobian;
                constraintLine];
            
            % Compute Kalman Gain
            K_k = P_k_k_1 * jacobian'...
                *((jacobian * P_k_k_1 * jacobian' + R) \ eye(C + 1));
            kalmanGainRecord(:,:, k) = K_k;
    
            % Obtain the innovation
            innovation = z_k - z_hat_k;
            innovationRecord(:, k) = innovation;
    
            % EKF's state update
            x_k_k = x_k_k_1 + K_k * innovation;
            x_k_k(WienerStatesSelection) = real(x_k_k(WienerStatesSelection));
     
            % EKF's covariance matrix update
            P_k_k = (eye(2*q + 1 + 4) - K_k*jacobian) * P_k_k_1;
            
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
        PStateRecord(:,:,k) = P_k_k;
        P_k_k_1 = F * P_k_k * F' + Q;
    
        errorStateRecord(:, k) = x_k_k(1:4);
        channelStateRecord(:, k) = x_k_k(5:end);
        if terminationCondition(x_LQG_k)
            fprintf("Terminated at step %d. Proceeding to save and proceed", k);
            earlyBreak = k;
            break;
        end
    
    end
    %% Plots
    lineWidth = 2;
    fontSize = 13;
    if earlyBreak ~= 0
        simulationSteps = earlyBreak;
    end
    epochVector = 1:simulationSteps;
    
    %Observe the STD of the innovation sequence time series\
    innovationStdRecord = std(innovationRecord(1:(end-1),1:simulationSteps),1,1);
    if  earlyBreak ~= 0
        innovationSum = Inf;
    else
        innovationSum = sum(innovationStdRecord);
    end

    figure(Name="STD of the innovations", NumberTitle="off");
    plot(epochVector, innovationStdRecord, 'LineWidth', lineWidth);
    ylabel("Standard deviation of the innovations");
    xlabel("Epochs (Simulation Steps)");
    
    % Observe the innovation sequence time series
    figure(Name="Middle tap of the innovation sequence", NumberTitle="off");
    hold on;
    plot(epochVector, real(innovationRecord(middleSample,1:simulationSteps)), 'LineWidth', lineWidth);
    plot(epochVector, imag(innovationRecord(middleSample,1:simulationSteps)), 'LineWidth', lineWidth);
    legend({"Real", "Imaginary"});
    ylabel("Innovation sequence of the middle tap");
    xlabel("Epochs (Simulation Steps)");
    hold off;
    
    figure(Name="Delay Estimation", NumberTitle="off");
    hold on;
    plot(epochVector, LQGStateRecord(1,1:simulationSteps), 'LineWidth', lineWidth);
    %plot(epochVector, LOSDelay(epochVector*4096), 'LineWidth', lineWidth);
    legend({"LQG's estimated delay", "True delay"});
    ylabel("Delay estimate");
    xlabel("Epochs (Simulation Steps)");
    hold off;
    
    figure(Name="Delay Error State", NumberTitle="off");
    hold on;
    plot(errorStateRecord(1,1:simulationSteps), 'LineWidth', lineWidth);
    plot(epochVector, zeros(1, simulationSteps), 'LineWidth', lineWidth);
    legend({"EKF's estimated delay error", "Zero line"});
    ylabel("Delay error estimate");
    xlabel("Epochs (Simulation Steps)");
    hold off;
    
    figure(Name="Phase Estimation", NumberTitle="off");
    hold on;
    plot(epochVector, LQGStateRecord(2,1:simulationSteps), 'LineWidth', lineWidth);
    %plot(epochVector, LOSPhase(epochVector*4096), 'LineWidth', lineWidth);
    legend({"LQG's estimated phase", "True Phase"});
    ylabel("Phase estimate");
    xlabel("Epochs (Simulation Steps)");
    hold off;
    
    figure(Name="Phase Error State", NumberTitle="off");
    hold on;
    plot(errorStateRecord(2,1:simulationSteps), 'LineWidth', lineWidth);
    plot(epochVector, zeros(1, simulationSteps), 'LineWidth', lineWidth);
    legend({"EKF's estimated phase error", "Zero line"});
    ylabel("Phase error estimate");
    xlabel("Epochs (Simulation Steps)");
    hold off;
    
    figure(Name="Doppler Estimation", NumberTitle="off");
    hold on;
    plot(epochVector, LQGStateRecord(3,1:simulationSteps)/(2*pi), 'LineWidth', lineWidth);
    %plot(epochVector, 2*pi*configuration.dopplerProfile(2) * ones(1,length(epochVector)), 'LineWidth', lineWidth);
    legend({"LQG's estimated Doppler frequency", "True Doppler frequency"});
    ylabel("Doppler estimate");
    xlabel("Epochs (Simulation Steps)");
    hold off;

    figure(Name="Trace of error matrix", NumberTitle="off");
    hold on;
    trace_vector = zeros(1, simulationSteps);
    for i = 1:simulationSteps
        trace_vector(1, i) = trace(PStateRecord(:,:,i));
    end
    plot(epochVector, trace_vector, 'LineWidth', lineWidth);
    %plot(epochVector, 2*pi*configuration.dopplerProfile(2) * ones(1,length(epochVector)), 'LineWidth', lineWidth);
    legend({"LQG's estimated Doppler frequency", "True Doppler frequency"});
    ylabel("Doppler estimate");
    xlabel("Epochs (Simulation Steps)");
    hold off;
    
    figure(Name="Doppler Error State", NumberTitle="off");
    hold on;
    plot(errorStateRecord(3,1:simulationSteps), 'LineWidth', lineWidth);
    plot(epochVector, zeros(1, simulationSteps), 'LineWidth', lineWidth);
    ylabel("Doppler error estimate");
    xlabel("Epochs (Simulation Steps)");
    legend({"EKF's estimated Doppler error", "Zero line"});
    hold off;
    
    figure(Name="Real Kalman Gain Elements", NumberTitle="off");
    hold on;
    for i = 1:size(kalmanGainRecord, 1)
        for j = 1:size(kalmanGainRecord, 2)
            plot(epochVector, squeeze(real(kalmanGainRecord(i, j,1:simulationSteps))), 'DisplayName', sprintf('K_{%d,%d}', i, j), 'LineWidth', lineWidth);
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
            plot(epochVector, squeeze(imag(kalmanGainRecord(i, j,1:simulationSteps))), 'DisplayName', sprintf('K_{%d,%d}', i, j), 'LineWidth', lineWidth);
        end
    end
    legend show;
    ylabel("Imaginary Kalman Gain Elements");
    xlabel("Epochs (Simulation Steps)");
    set(gca, "FontSize", fontSize);
    hold off;
    
    figure(Name="Channel Weights Over Time", NumberTitle="off");
    hold on;
    for i = 1:3
        plot(abs(channelStateRecord(i,1:simulationSteps)), 'LineWidth', lineWidth);
    end
    yyaxis right
    plot(constraintRecord, 'LineWidth', lineWidth);
    ylabel("Constraint");
    xlabel("Epochs (Simulation Steps)");
    hold off;
    
    figure(Name="", NumberTitle="off");
    hold on;
    for i = 1:q
        plot(abs(channelStateRecord(i+q+1,1:simulationSteps)), 'LineWidth', lineWidth);
    end
    hold off;

    figure(Name="main tap", NumberTitle="off");
    hold on; 
    plot(real(channelStateRecord(q+1,1:simulationSteps)), 'LineWidth', lineWidth);
    hold off;
    %% Preparing folder
    % Right now we save the results of the run on its own folder
    %% 1. Global Folder Management
    % Check if the global directory exists; if not, create it
    if ~exist(parentFolderName, 'dir')
        mkdir(parentFolderName);
        fprintf('Initial run detected. Created global folder: %s\n', parentFolderName);
    end
    
    % Add the folder to the path so MATLAB can find files within it easily
    addpath(genpath(parentFolderName));

    %% 2. Individual Run Folder
    % Create a subfolder inside the parent for this specific run
    timestamp = int2str(currentRun);
    runFolderName = ['SimRun_', timestamp];
    
    % Construct the full path (e.g., "MySimulations/SimRun_20260428_210000")
    runPath = fullfile(parentFolderName, runFolderName);
    mkdir(runPath);
    % Save the configuration used for this specific run
    save(fullfile(runPath, 'config_used.mat'), 'configuration');
    save(fullfile(runPath, 'initiation_used.mat'), 'initialCondition');

    %% Automated Image Saving
    % Find all open figure handles
    figHandles = findobj('Type', 'figure');
    
    for i = 1:length(figHandles)
        % Get the name of the figure to use as a filename
        figName = get(figHandles(i), 'Name');
        if isempty(figName)
            figName = sprintf('Figure_%d', figHandles(i).Number);
        end
        
        % Clean the filename (remove special characters/spaces)
        cleanName = regexprep(figName, '[^a-zA-Z0-9_]', '_');
        
        % Save as PNG (better for quick viewing)
        saveas(figHandles(i), fullfile(runPath, [cleanName, '.png']));
        
        % Optional: Save as FIG (better for post-analysis)
        savefig(figHandles(i), fullfile(runPath, [cleanName, '.fig']));
    end
    
    fprintf('Simulation %i completed', currentRun);