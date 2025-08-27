function [errorCovarianceAPriori] = propagate_covariance(currentErrorCovariance, modelTransitionMatrix, stateNoiseCovariance)
%   PROPAGATE_COVARIANCE 
%   calculates the covariance a priori given the current covariance, the state model
%   transition matrix and the state noise covariance

errorCovarianceAPriori = ...
modelTransitionMatrix * currentErrorCovariance * conj(modelTransitionMatrix)' ...
+ stateNoiseCovariance;
end

