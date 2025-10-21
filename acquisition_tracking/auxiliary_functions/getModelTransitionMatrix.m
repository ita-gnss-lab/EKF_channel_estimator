function [carrierStateTransitionMatrix, channelStateTransitionMatrix] = getModelTransitionMatrix(epoch, carrierFrequency, numberOfTaps)
%GETMODELTRANSITIONMATRIX Summary of this function goes here
%   Detailed explanation goes here
beta = 1 / (2 * pi * carrierFrequency);
T = epoch;

carrierStateTransitionMatrix = [
1 0 beta*T 0.5*beta*T^2;
0 1 T      0.5*T^2;
0 0 1      T;
0 0 0      1
];
channelStateTransitionMatrix = eye(numberOfTaps + 1);
end
