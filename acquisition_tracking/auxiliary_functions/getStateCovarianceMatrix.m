function [stateCovarianceMatrix] = getStateCovarianceMatrix(varianceSquared, epoch, carrierFrequency, numberOfTaps)
%GETCOVARIANCEMATRIX Summary of this function goes here
%   Detailed explanation goes here
beta = 1 / (2 * pi * carrierFrequency);
T = epoch;

zeroOrder = zeros(4);
zeroOrder(1,1) = T;
firstOrder = zeros(4);
firstOrder(2,2) = T;
secondOrder = [
(T^3 * beta^2 / 3) (T^3 * beta / 3) (T^2 * beta / 2) 0;
(T^3 * beta / 3)   (T^3 / 3)        (T^2 / 2)        0;
(T^2 * beta / 2)   (T^2 / 2)         T               0;
zeros(1, 4)
];
thirdOrder = [
(T^5 * beta^2 / 20) (T^5 * beta / 20) (T^4 * beta / 8) (T^3 * beta / 6);
(T^5 * beta / 20)   (T^5 / 20)        (T^4 / 8)        (T^3 / 6);
(T^4 * beta / 8)    (T^4 / 8)         (T^3 / 3)        (T^2 / 2);
(T^3 * beta / 6)    (T^3 / 6)         (T^2 / 2)         T
];

errorStateCovariance = varianceSquared(1) * zeroOrder ...
    + varianceSquared(2) * firstOrder ...
    + varianceSquared(3) * secondOrder ...
    + varianceSquared(4) * thirdOrder;

channelStateCovariance = varianceSquared(5) * eye(numberOfTaps + 1);

stateCovarianceMatrix = blkdiag(errorStateCovariance, channelStateCovariance);
end

