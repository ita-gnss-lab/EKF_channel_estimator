clear all;
close all;

% Parameters
    % Simulation
    K = 500;

    % System
    carrierFrequency = 1575.42e6;
    beta = 1/(2*pi*carrierFrequency);
    sequencePeriod = 1e-3;
    % Shouldn't the sampling frequency be a multiple of 1.023e6, due to the chip rate?
    samplingFrequency = 4e6;
    samplingPeriod = 1/samplingFrequency;
    % Maybe we should rename this to chipRate. Symbol is often referred to 
    % as the navigation data bit, which should not be adopted here.
    symbolFrequency = 1.023e6;
    % chipPeriod?
    symbolPeriod = 1/symbolFrequency;
    numberOfPeriods = 1;
    samplesPerPeriod = samplingFrequency * sequencePeriod;
    samplesTotal = numberOfPeriods * samplesPerPeriod;
    inputSampleIndex = 100 * samplesTotal;
    
    % Channel
    numberOfTaps = 10;
    
    % Variances
    varianceDoppDrift = 0.02;%3.141e-4;
    varianceClockN1 = 5.141e-5;%3.141e-4;
    varianceClockN2 = 5.141e-5;%3.141e-4;
    varianceClockN3 = 5.141e-5;%3.141e-4;
    varianceIonosphere = 5.141e-5;%3.141e-4;

    % Math
    Circulant = Circulant_Toeplitz_Matrices(samplesTotal, numberOfTaps);

% Model

channel = ones(1, numberOfTaps + 1);
modelTransitionMatrix = [
1   0   beta*sequencePeriod   beta*sequencePeriod^2/2   zeros(1, numberOfTaps + 1);
0   1   sequencePeriod        sequencePeriod^2/2        zeros(1, numberOfTaps + 1);
0   0   1                     sequencePeriod            zeros(1, numberOfTaps + 1);
0   0   0                     1                         zeros(1, numberOfTaps + 1);
zeros(numberOfTaps + 1, 4)                                  diag(channel)
];

clockCovariance = ...
varianceClockN1 * ...
[
beta^2*sequencePeriod beta*sequencePeriod 0 0
beta*sequencePeriod   sequencePeriod      0 0 
zeros(2, 4)
] + ...
varianceClockN2 * ...
[
beta^2*sequencePeriod^3/3 beta*sequencePeriod^3/3 beta*sequencePeriod^2/2 0;
beta*sequencePeriod^3/3   sequencePeriod^3/3      sequencePeriod^2/2      0;
beta*sequencePeriod^2/2   sequencePeriod^2/2      sequencePeriod          0;
zeros(1,4)
] + ...
varianceClockN3 * ...
[
beta^2*sequencePeriod^5/20 beta*sequencePeriod^5/20 beta*sequencePeriod^4/8 beta*sequencePeriod^4/6
beta*sequencePeriod^5/20   sequencePeriod^5/20      sequencePeriod^4/8      sequencePeriod^3/6
beta*sequencePeriod^4/8    sequencePeriod^4/8       sequencePeriod^3/3      sequencePeriod^2/2
beta*sequencePeriod^3/6    sequencePeriod^3/6       sequencePeriod^2/2      sequencePeriod
];

angDoppDriftCovariance = ...
varianceDoppDrift * ...
[
beta^2*sequencePeriod^5/20 beta*sequencePeriod^5/20 beta*sequencePeriod^4/8 beta*sequencePeriod^4/6
beta*sequencePeriod^5/20   sequencePeriod^5/20      sequencePeriod^4/8      sequencePeriod^3/6
beta*sequencePeriod^4/8    sequencePeriod^4/8       sequencePeriod^3/3      sequencePeriod^2/2
beta*sequencePeriod^3/6    sequencePeriod^3/6       sequencePeriod^2/2      sequencePeriod
];

ionosphereCovariance = ...
varianceIonosphere * ...
[
sequencePeriod 0 0 0;
0              0 0 0;
0              0 0 0;
0              0 0 0
];

firstStateNoiseCovariance = clockCovariance + angDoppDriftCovariance + ionosphereCovariance;

stateNoiseCovariance = ...
    [
    firstStateNoiseCovariance zeros(4, numberOfTaps + 1);
    zeros(numberOfTaps + 1, 4) diag(0.01*ones(1, numberOfTaps+1))
    ];

% State-Vector

stateDelay          = [];
stateDoppPhase      = [];
stateDopp           = [];
stateDoppDrift      = [];
initChannelWeights  = [0.2 zeros(1, numberOfTaps)];
quadrature          = [];
inPhase             = [];
stateError          = [];
changeSignal        = [zeros(15, 1)];
changeSignalIntegration = [zeros(15, 1)];

% Covariances
doppPhaseVariance = pi^2/3;                    % Phase Variance
% Doppler Variance
% NOTE: I think this is too large. Acquisition often give a 
% frequency estimate with +-25 Hz precision. Assuming the Doppler frequency
% estimate as a random variable with uniform probability density function 
% with bounds [a, b], we have that Var[f_D] = (b-a)^2 / 12. Thus, we would
% have that doppVariance = 50^2 / 12.
% NOTE: From my experiments with Kaman filter based PLLs, i think that very
% large initial Doppler frequency shift variances could make the filter to
% never reach convergence.
doppVariance = 250^2*4*pi^2/3; 
doppDriftVariance = 0.02;                 % Doppler Drift Variance
initDelayDoppCovariance = [
0.5*symbolPeriod 0 0 0;
0 doppPhaseVariance 0 0; 
 0 0 doppVariance 0;
0 0 0 doppDriftVariance];
initChanCovariance = eye(numberOfTaps + 1)*0.1;

initCovariance = [
initDelayDoppCovariance zeros(4, numberOfTaps + 1)    
zeros(numberOfTaps + 1, 4) initChanCovariance
];


% Measurement Sensibility Bases

[realDelaySensibilityBase, imagDelaySensibilityBase] = delayLinearization(numberOfTaps, symbolPeriod);
realDelaySensibilityBase = realDelaySensibilityBase/samplesTotal;
imagDelaySensibilityBase = imagDelaySensibilityBase/samplesTotal;

[realDoppDriftSensibilityBase, imagDoppDriftSensibilityBase] = doppDriftLinearization(numberOfTaps, symbolPeriod);
realDoppDriftSensibilityBase = realDoppDriftSensibilityBase/samplesTotal;
imagDoppDriftSensibilityBase = imagDoppDriftSensibilityBase/samplesTotal;

[realAmbiguityVectorBase, imagAmbiguityVectorBase] = ambiguityVector(numberOfTaps, symbolPeriod);
realAmbiguityVectorBase = realAmbiguityVectorBase/samplesTotal;
imagAmbiguityVectorBase = imagAmbiguityVectorBase/samplesTotal;
%Signal File Name

fileName = 'signal_source_L1E1_GNSSR_2.dat';


%% Tracking and Acquisition

for k = 1 : K

    [signal, outputSampleIndex, File_Ended] = read_gr_complex_binary ( ...
        fileName, ...
        inputSampleIndex, ...
        samplesTotal);
    
    if File_Ended == true
        break;
    end
    disp(inputSampleIndex);

    inputSampleIndex = outputSampleIndex;

    % Time Vector
    
    time = 0 : samplingPeriod : (samplesTotal - 1)*samplingPeriod;

    if k <= 1
    
        % Signal Acquisition
        if exist('ACQ_DATA_L1E1_GNSSR_2.mat') ~= 2
            Number_Samples_per_Period = samplingFrequency*sequencePeriod;
            [ACQ_DATA, Doppler_Vector, Threshold, Satellites] = ca_acquisition( ...
                time(1 : Number_Samples_per_Period), ...
                signal(1 : Number_Samples_per_Period), ...
                samplingFrequency, ...
                symbolFrequency, ...
                1);

            save ACQ_DATA_L1E1_GNSSR_2.mat ACQ_DATA Doppler_Vector Threshold Satellites;
        
        else

            load ACQ_DATA_L1E1_GNSSR_2.mat;
            %Satellites = sats_found;
            %Doppler_Vector = doppler_bin_vec;

        end
        
        disp(Satellites);

        % Delay and Doppler Acquisition

        for Satellite = 1 %: length(Satellites)

            stateDelay                           = [stateDelay        ACQ_DATA(Satellites(Satellite)).max_index(2)*samplingPeriod];
            stateDoppPhase                       = [stateDoppPhase    0];
            stateDopp                       = [stateDopp    Doppler_Vector(ACQ_DATA(Satellites(Satellite)).max_index(1))];
            stateDoppDrift                    = [stateDoppDrift 0];
            channelWeights(:, 1, Satellite)      = [0.2 zeros(1, numberOfTaps)];
            stateError(:, 1, Satellite)          = zeros(2*numberOfTaps + 1, 1);
            errorCovariance(1, Satellite, :, :)  = initCovariance;
            quadrature = [quadrature 0];
            inPhase = [inPhase 0];
            

        end
    
    else

        % Tracking

        for Sattellite = 1 %: length(Satellites)

            % Kalman Propagation Step 
            % Current state
            stateVector = [
                            stateDelay(k-1, Sattellite);
                            stateDoppPhase(k-1, Sattellite);
                            2*pi*stateDopp(k-1, Sattellite);
                            2*pi*stateDoppDrift(k-1, Sattellite);
                            channelWeights(:, k-1, Sattellite)
                           ];

            currentErrorCovariance(:,:) = errorCovariance(k-1, Sattellite, :, :);

            % Propagation Step

            stateAPriori = modelTransitionMatrix * stateVector;
            % NOTE: I think that it is not necessary to use
            % `conj(modelTransitionMatrix)` here, given that
            % `modelTransitionMatrix` would already yield an Hermitian
            % matrix.
            % NOTE: Note also that modelTransitionMatrix is completly real,
            % so it would be sufficient to use `modelTransitionMatrix.'`
            % here.
            % NOTE: In general, the Kalman filter is defined for real
            % numbers. Andreas Iliopoulos, however, adopt a complex Kalman
            % filter for tracking the channel coefficients for each
            % multipath signal as complex numbers.
            errorCovarianceAPriori = ...
            modelTransitionMatrix * currentErrorCovariance * conj(modelTransitionMatrix)' ...
            + stateNoiseCovariance;

            % Doppler Correction 
            stateDoppPhaseAPriori = stateAPriori(2);
            stateDoppAPriori = stateAPriori(3);
            dopplerEvolution = exp(1i*(stateDoppPhaseAPriori + time*stateDoppAPriori));
            processedSignal = signal; % I did not correct for Doppler
            % dopp = exp(1i*(Phase_Satellites_Matrix(Interval,Sat)+2*pi*Doppler_Satellites_Matrix(Interval,Sat).*Time))
            % Signal_Corrected = conj(conj(Measured_Signal).*dopp)';

            % Correlation
            stateDelayAPriori = real(stateAPriori(1));
            Multi_Correlator = [];
            disp(stateDelayAPriori)
            for aux = -numberOfTaps : 1 : numberOfTaps
                Multi_Correlator = [Multi_Correlator; 
                reference_signal(Satellites(Satellite), ...
                stateDelayAPriori + aux*symbolPeriod/numberOfTaps, ...
                symbolFrequency, ...
                samplingFrequency, ...
                numberOfPeriods).*dopplerEvolution];
            end
            measurementNoiseCovariance = Multi_Correlator*Multi_Correlator';
            % Real Measure
            Z = conj(processedSignal) * conj(Multi_Correlator)' / samplesTotal;
            % stem(abs(Z));
            % pause(0.1);
            
            % Estimated Measure
            channelWeightsAPriori = stateAPriori(5:end);
            channelMatrix = zeros(samplesTotal);
            for L = 1:numberOfTaps
                channelMatrix = channelMatrix + channelWeightsAPriori(L)*Circulant(:, :, L);
            end
            
            Ref = reference_signal(Satellites(Satellite), ...
                stateDelayAPriori, ...
                symbolFrequency, ...
                samplingFrequency, ...
                numberOfPeriods)';
            
            % Estimated_Measurement = Multi_correlator*conj(Estimated_Channel_Matrix*Estimated_Sequence)/ Number_of_Samples;
            signalAPriori = channelMatrix * (Ref).* conj(dopplerEvolution)';
            ZAPriori = signalAPriori'*conj(Multi_Correlator)'/samplesTotal;

            % Error Signal
            errorSignal = Z - ZAPriori;
            stateError(:, k, Sattellite) = errorSignal;

            %% Measurement Sensibility Matrix
            sensibilityMatrix = [];

            % Delay
            delaySensibility = zeros(2*numberOfTaps + 1, 1);
            auxiliaryValue = 1;%exp(-1j*stateDoppPhaseAPriori);
            for L = 1 : numberOfTaps
               Linearization = auxiliaryValue*conj(channelWeightsAPriori(L))*(realDelaySensibilityBase(:, L) + 1j*imagDelaySensibilityBase(:, L)); 
               delaySensibility  = delaySensibility + Linearization;
            end
            sensibilityMatrix = [sensibilityMatrix delaySensibility];
            
            % Doppler Phase
            doppPhaseSensibility = zeros(2*numberOfTaps + 1, 1);
            auxiliaryValue = -1j; %*exp(-1j*stateDoppPhaseAPriori);
            for L = 1 : numberOfTaps
               Linearization = auxiliaryValue*conj(channelWeightsAPriori(L))*(realAmbiguityVectorBase(:, L) + 1j*imagAmbiguityVectorBase(:, L)); 
               doppPhaseSensibility  = doppPhaseSensibility + Linearization;
            end
            sensibilityMatrix = [sensibilityMatrix doppPhaseSensibility];
            
            % Doppler Drift
            doppDriftSensibility = zeros(2*numberOfTaps + 1, 1);
            auxiliaryValue = 1; %exp(-1j*stateDoppPhaseAPriori);
            for L = 1 : numberOfTaps
               Linearization = auxiliaryValue*conj(channelWeightsAPriori(L))*(realDoppDriftSensibilityBase(:, L) + 1j*imagDoppDriftSensibilityBase(:, L)); 
               doppDriftSensibility  = doppDriftSensibility + Linearization;
            end
            sensibilityMatrix = [sensibilityMatrix doppDriftSensibility];

            % Angular Doppler 
            angDoppSensibility = zeros(2*numberOfTaps + 1, 1);
            sensibilityMatrix = [sensibilityMatrix angDoppSensibility];

            % Channel Coefficients
            channelSensibility = zeros(2*numberOfTaps + 1, numberOfTaps + 1);
            auxiliaryValue = 1; %exp(-1j*stateDoppPhaseAPriori);
            for L = 1 : numberOfTaps+1
               Linearization = auxiliaryValue*(realAmbiguityVectorBase(:, L) + 1j*imagAmbiguityVectorBase(:, L)); 
               channelSensibility(:,L)  = Linearization;
            end
            sensibilityMatrix = [sensibilityMatrix channelSensibility];

            %% Kalman Estimation

            % Kalman gain
            inverse = inv(sensibilityMatrix*errorCovarianceAPriori*conj(sensibilityMatrix)' + measurementNoiseCovariance);
            kalmanGain = errorCovarianceAPriori*conj(sensibilityMatrix)'*(inverse);

            % State Estimate Update
            changeImpulse = kalmanGain*conj(errorSignal)';
            changeSignal = [changeSignal changeImpulse];
            changeSignalIntegration = [changeSignalIntegration (changeSignalIntegration(end) + changeImpulse)];
            stateAPosteriori = stateAPriori + changeImpulse;

            % Error Covariance Update
            errorCovarianceAPosteriori = (eye(numberOfTaps + 1 + 4) - kalmanGain*sensibilityMatrix)*errorCovarianceAPriori;

            % In-Phase and Quadrature
            centralPrompt = reference_signal(Satellites(Satellite), ...
                real(stateAPosteriori(1)), ...
                symbolFrequency, ...
                samplingFrequency, ...
                numberOfPeriods);
            quadrature(k-1, Sattellite)=imag(conj(processedSignal) * centralPrompt' / samplesTotal);
            inPhase(k-1, Sattellite)=real(conj(processedSignal) * centralPrompt' / samplesTotal);
            

            % Make new current state
            stateDelay(k, Sattellite) = real(stateAPosteriori(1));
            stateDoppPhase(k, Sattellite) = real(stateAPosteriori(2));
            stateDopp(k, Sattellite) = real(stateAPosteriori(3)) / (2*pi);
            stateDoppDrift(k, Sattellite) = real(stateAPosteriori(4)) / (2*pi);
            channelWeights(:, k, Sattellite) = stateAPosteriori(5:end);
            errorCovariance(k, Sattellite, :, :) = errorCovarianceAPosteriori;
        end      
    end 
end

figure('Name','Delay Error','NumberTitle','off');
plot(real(stateError(1, :, 1)));grid on;

figure('Name','Delay','NumberTitle','off');
plot(stateDelay(:, 1));grid on;



% figure('Name','Phase error Mine','NumberTitle','off');
% plot(Phase_error(:, 1));grid on;

figure('Name','Doppler Phase','NumberTitle','off');
plot(stateDoppPhase(:, 1));grid on;


% figure('Name','Doppler error Mine','NumberTitle','off');
% plot(Doppler_error(:, 1));grid on;

figure('Name','Doppler','NumberTitle','off');
plot(stateDopp(:, 1));grid on;


figure('Name','Doppler Drift','NumberTitle','off');
plot(stateDoppDrift(:, 1));grid on;

figure('Name','Main reflection','NumberTitle','off');
plot(abs(channelWeights(1, :, 1)));grid on;

 figure('Name','In-Phase Mine','NumberTitle','off');
 plot(inPhase(:, 1));grid on;

 figure('Name','Quadrature Mine','NumberTitle','off');
 plot(quadrature(:, 1));grid on;
% 
% figure('Name', 'C/N0 Mine' , 'NumberTitle','off');
% plot(ctnd(:, 1));grid on;

