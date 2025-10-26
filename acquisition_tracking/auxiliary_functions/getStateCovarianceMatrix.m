function [stateCovarianceMatrix] = getStateCovarianceMatrix(sigma2Vec, epoch, beta, numberOfTaps)
%GETCOVARIANCEMATRIX Summary of this function goes here
%   Detailed explanation goes here
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

errorStateCovariance = sigma2Vec(1) * zeroOrder ...
    + sigma2Vec(2) * firstOrder ...
    + sigma2Vec(3) * secondOrder ...
    + sigma2Vec(4) * thirdOrder;

channelStateCovariance = blkdiag(sigma2Vec(5), sigma2Vec(6) * eye(numberOfTaps));

stateCovarianceMatrix = blkdiag(errorStateCovariance, channelStateCovariance);
end
