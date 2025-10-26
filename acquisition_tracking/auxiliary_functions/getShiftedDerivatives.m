function shiftedDerivatives = getShiftedDerivatives(delayError, q, configuration)
%GETSHIFTEDDERIVATIVES Builds dPhi_pp/d(eps_tau) over taps/correlators.
%
%   delayError : scalar eps_tau
%   q          : correlator half-span (m in {q,...,-q})
%   configuration : receiver configuration w/ sampling & chipping freq
%
% Returns a (2q+1)x(q+1) matrix whose (m,l) entry corresponds to the derivative
% of the autocorrelation evaluated at Delta_{m,l} = eps_tau + (l + m)Ts. The row
% order (top-to-bottom) follows m = q, q-1, ..., -q so it matches
% buildCorrelatorBank().

    Ts = 1 / configuration.samplingFrequency;
    m = (q:-1:-q).';      % (2q+1)x1
    l = 0:q;              % 1x(q+1)

    deltaMatrix = delayError + (l + m) * Ts;

    shiftedDerivatives = zeros(size(deltaMatrix));
    for i = 1:size(deltaMatrix, 1)
        for j = 1:size(deltaMatrix, 2)
            shiftedDerivatives(i, j) = autocorrelationDerivative( ...
                deltaMatrix(i, j), ...
                configuration.chippingFrequency, ...
                configuration.samplingFrequency);
        end
    end
end
