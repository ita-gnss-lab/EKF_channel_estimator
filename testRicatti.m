%% Parameters
simulationSteps = 500;
q = 7;
C = 2*q + 1;
middleSample = q + 1;
epoch = configuration.totalChips / configuration.chippingFrequency;

% NOTE: These were not being used
WienerStatesSelection = 1:4;
% channelWeights = 5 + 0:numberOfTaps;

%% Covariances 
% Convert CN0 from dB-Hz to linear scale
carrierToNoiseRatioLinear = 10^(configuration.carrierToNoiseDensityRatio / 10);
% Compute the noise variance
thermalNoiseVarianceSquared = configuration.samplingFrequency / carrierToNoiseRatioLinear;

sigma2WVec = [1e-2 0 0 0 0];
Q = getStateCovarianceMatrix(...
    sigma2WVec, ...
    epoch, ...
    configuration.carrierFrequency, ...
    q);

%% State History Vectors
LQGStateRecord = zeros(4, simulationSteps);
errorStateRecord = zeros(4, simulationSteps);
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
T_e =  blkdiag(beta, 1, 1/epoch, 2/epoch^2);
T_u =  blkdiag(beta, 1, 1/epoch, 2/epoch^2);

%% Coupling Matrix for Control Signal
B_LQG = eye(4);

%% IDARE Solution
maxIter = 100;
traceL = zeros(maxIter,1);
relation = (1:maxIter) * 0.000001;

for i = (1:maxIter)
    [~, L, ~] = idare(F_W, ...
    B_LQG, ...
    T_e * relation(i), ...
    T_u, ...
    [], []);
    traceL(i) = trace(L);
end

figure;
plot(traceL);


